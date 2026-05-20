"""Single mutation pipeline — the core REFOLD inference loop.

Steps:
1. Fetch AlphaFold structure + parse
2. Apply mutation + predict ΔΔG (ESM-1v ensemble)
3. Rescue amenability classification
4. ANM ensemble + fpocket transient pocket detection
5. De novo molecule generation conditioned on transient pockets
6. Filter, score, rank top candidates
"""

from __future__ import annotations

import argparse
import logging
import time
from functools import cached_property
from pathlib import Path
from typing import Optional

import numpy as np

from refold.constants import (
    CHECKPOINT_DIR, RESULTS_DIR,
    DDG_DESTABILIZING_THRESHOLD, DDG_SEVERELY_UNSTABLE,
    RESCUE_PROBABILITY_THRESHOLD, HIGH_CONFIDENCE_THRESHOLD,
    N_CONFORMATIONS, N_MOLECULES_PER_POCKET,
    ESM2_EMBED_DIM, ESM1V_MAX_SEQ_LEN,
    GNN_NODE_DIM, GNN_EDGE_DIM, THERMO_FEAT_DIM, EVO_FEAT_DIM,
    CONTACT_CA_CUTOFF,
)
from refold.types import (
    Mutation, MutationClass, RescueAmenability,
    ProteinStructure, REFOLDResult, GeneratedMolecule,
)
from refold.exceptions import StructureNotFoundError

logger = logging.getLogger(__name__)


class REFOLDPipeline:
    """End-to-end REFOLD pipeline for a single missense mutation."""

    def __init__(
        self,
        checkpoint_dir: Optional[Path] = None,
        device: "Optional[torch.device]" = None,
        n_conformations: int = N_CONFORMATIONS,
        n_molecules_per_pocket: int = N_MOLECULES_PER_POCKET,
        n_diffusion_steps: Optional[int] = None,
    ):
        self._checkpoint_dir = checkpoint_dir or CHECKPOINT_DIR
        self._device_arg = device
        self.n_conformations = n_conformations
        self.n_molecules_per_pocket = n_molecules_per_pocket
        self.n_diffusion_steps = n_diffusion_steps

    @cached_property
    def device(self) -> "torch.device":
        if self._device_arg is not None:
            return self._device_arg
        from refold.utils.device import get_device
        return get_device()

    @cached_property
    def esm1v_predictor(self) -> "ESM1vStabilityPredictor":
        from refold.structure.mutator import ESM1vStabilityPredictor
        return ESM1vStabilityPredictor(device=self.device)

    @cached_property
    def structure_predictor(self) -> "StructureBasedDDGPredictor":
        from refold.structure.mutator import StructureBasedDDGPredictor
        ckpt = str(self._checkpoint_dir / "rescue_classifier" / "best.pt")
        return StructureBasedDDGPredictor(checkpoint_path=ckpt, device=self.device)

    @cached_property
    def rescue_classifier(self) -> "Optional[RescueClassifier]":
        from refold.models.rescue_classifier.model import RescueClassifier
        import torch
        ckpt_path = self._checkpoint_dir / "rescue_classifier" / "best.pt"
        if not ckpt_path.exists():
            logger.warning(
                f"Rescue classifier checkpoint not found at {ckpt_path} — "
                "using heuristic classification"
            )
            return None
        try:
            model = RescueClassifier()
            ckpt = torch.load(ckpt_path, map_location=self.device)
            model.load_state_dict(ckpt.get("model_state_dict", ckpt), strict=False)
            model = model.to(self.device)
            model.eval()
            return model
        except Exception as e:
            logger.warning(f"Failed to load rescue classifier: {e}")
            return None

    @cached_property
    def molecule_generator(self) -> "Optional[REFOLDDiffusionModel]":
        from refold.models.molecule_generator.diffusion import REFOLDDiffusionModel
        import torch
        ckpt_path = self._checkpoint_dir / "molecule_generator" / "best.pt"
        T = self.n_diffusion_steps or 1000
        model = REFOLDDiffusionModel(T=T)
        if not ckpt_path.exists():
            logger.warning(
                f"Molecule generator checkpoint not found at {ckpt_path} — "
                "using untrained model (research only)"
            )
            return model.to(self.device)
        try:
            ckpt = torch.load(ckpt_path, map_location=self.device)
            model.load_state_dict(ckpt.get("model_state_dict", ckpt), strict=False)
            return model.to(self.device)
        except Exception as e:
            logger.warning(f"Failed to load molecule generator: {e}")
            return model.to(self.device)

    def run(
        self,
        mutation: Mutation,
        output_dir: Optional[Path] = None,
    ) -> REFOLDResult:
        """Execute the full REFOLD pipeline for a single missense mutation."""
        t_start = time.perf_counter()

        if output_dir is None:
            output_dir = RESULTS_DIR / "molecule_candidates" / mutation.uniprot_id

        error_message: Optional[str] = None
        mutation_class = MutationClass.UNKNOWN
        rescue_amenability = RescueAmenability.LOW
        rescue_prob = 0.0
        ddg = 0.0
        pockets: list = []
        all_molecules: list[GeneratedMolecule] = []
        top_candidates: list[GeneratedMolecule] = []

        try:
            # ── Step 1: Fetch structure ─────────────────────────────────────
            logger.info(f"[{mutation}] Step 1: Fetching AlphaFold structure")
            from refold.data.alphafold_fetcher import fetch_structure
            from refold.structure.parser import parse_pdb_to_structure

            pdb_path = fetch_structure(mutation.uniprot_id)
            if pdb_path is None:
                raise StructureNotFoundError(mutation.uniprot_id)

            structure = parse_pdb_to_structure(pdb_path, mutation.uniprot_id)

            # ── Step 2: Apply mutation + ΔΔG ──────────────────────────────
            logger.info(f"[{mutation}] Step 2: Predicting ΔΔG")
            from refold.structure.mutator import apply_mutation, predict_ddg_ensemble

            mutant_struct = apply_mutation(structure, mutation)
            ddg, ddg_esm1v, ddg_gnn = predict_ddg_ensemble(
                structure.sequence, structure, mutation,
                self.esm1v_predictor, self.structure_predictor,
            )
            mutant_struct.ddg_predicted = ddg
            logger.info(f"[{mutation}] ΔΔG = {ddg:.2f} kcal/mol")

            # ── Step 3: Rescue amenability classification ──────────────────
            logger.info(f"[{mutation}] Step 3: Rescue amenability classification")
            esm2_emb = self._get_esm2_embedding(structure.sequence)
            rescue_prob = self._classify_with_model(structure, mutation, ddg, esm2_emb)
            mutation_class, rescue_amenability = self._classify_mutation(ddg, rescue_prob)

            logger.info(
                f"[{mutation}] Class={mutation_class.value} "
                f"Amenability={rescue_amenability.value} "
                f"prob={rescue_prob:.3f}"
            )

            if rescue_amenability == RescueAmenability.UNRESCUABLE:
                logger.info(f"[{mutation}] Mutation unrescuable — skipping pocket detection")
            else:
                # ── Step 4: Conformational ensemble + pocket detection ─────
                logger.info(f"[{mutation}] Step 4: Generating conformational ensemble")
                from refold.structure.ensemble import (
                    generate_conformational_ensemble, generate_mutant_ensemble,
                )
                from refold.structure.pocket_finder import detect_transient_pockets

                wt_conformations = generate_conformational_ensemble(
                    structure, n_conformations=self.n_conformations
                )
                mutant_conformations = generate_mutant_ensemble(
                    mutant_struct, n_conformations=self.n_conformations
                )
                pockets = detect_transient_pockets(
                    structure, wt_conformations, mutant_conformations,
                    full_pdb_path=pdb_path,
                )
                n_transient = sum(1 for p in pockets if p.is_transient)
                logger.info(f"[{mutation}] {len(pockets)} pockets ({n_transient} transient)")

                if pockets and self.n_molecules_per_pocket > 0:
                    # ── Step 5: Molecule generation ────────────────────────
                    logger.info(f"[{mutation}] Step 5: Generating molecules")
                    gen = self.molecule_generator
                    if gen is not None:
                        from refold.models.molecule_generator.decoder import decode_molecule
                        target_pockets = sorted(
                            pockets,
                            key=lambda p: (p.is_transient, p.druggability_score),
                            reverse=True,
                        )[:3]

                        for pocket in target_pockets:
                            raw_samples = gen.sample(
                                pocket,
                                n_molecules=self.n_molecules_per_pocket,
                                device=self.device,
                            )
                            for coords, atom_types in raw_samples:
                                smiles = decode_molecule(coords, atom_types)
                                if smiles:
                                    mol = GeneratedMolecule(
                                        smiles=smiles,
                                        pocket_id=pocket.pocket_id,
                                    )
                                    all_molecules.append(mol)

                        # ── Step 6: Filter, score, rank ────────────────────
                        logger.info(
                            f"[{mutation}] Step 6: Filtering "
                            f"{len(all_molecules)} molecules"
                        )
                        from refold.scoring.filters import filter_molecule_library
                        from refold.scoring.binding_affinity import (
                            estimate_binding_affinity_heuristic,
                        )
                        from refold.scoring.rescue_probability import (
                            compute_final_rescue_probability,
                            rank_molecules_by_rescue_probability,
                        )

                        pocket_map = {m.smiles: p for p in target_pockets
                                      for m in all_molecules if m.pocket_id == p.pocket_id}

                        passing = filter_molecule_library(all_molecules)

                        for mol in passing:
                            p = pocket_map.get(mol.smiles, target_pockets[0])
                            mol.predicted_affinity_kcal = (
                                estimate_binding_affinity_heuristic(mol, p)
                            )
                            mol.rescue_probability = compute_final_rescue_probability(
                                mol, p, rescue_prob, ddg
                            )

                        top_candidates = rank_molecules_by_rescue_probability(passing)[:10]

                        logger.info(
                            f"[{mutation}] {len(all_molecules)} generated → "
                            f"{len(passing)} pass filters → "
                            f"top: {top_candidates[0].smiles if top_candidates else 'none'}"
                        )

        except StructureNotFoundError as e:
            error_message = str(e)
            logger.error(f"[{mutation}] {e}")
        except Exception as e:
            error_message = str(e)
            logger.exception(f"[{mutation}] Pipeline error: {e}")

        runtime = time.perf_counter() - t_start

        result = REFOLDResult(
            mutation=mutation,
            mutation_class=mutation_class,
            rescue_amenability=rescue_amenability,
            rescue_amenability_prob=rescue_prob,
            rescue_probability=rescue_prob,
            ddg_predicted=ddg,
            n_pockets_detected=len(pockets),
            n_molecules_generated=len(all_molecules),
            n_molecules_passing_filters=len([m for m in all_molecules if m.passes_all_filters]),
            top_candidates=top_candidates,
            pockets=pockets,
            all_molecules=all_molecules,
            error_message=error_message,
            runtime_seconds=runtime,
        )

        self._save_result(result, output_dir)
        return result

    def _classify_mutation(
        self,
        ddg: float,
        rescue_prob: float,
    ) -> tuple[MutationClass, RescueAmenability]:
        """Classify mutation class and rescue amenability."""
        if ddg > DDG_SEVERELY_UNSTABLE:
            return MutationClass.MISFOLDING, RescueAmenability.UNRESCUABLE
        if ddg < DDG_DESTABILIZING_THRESHOLD:
            return MutationClass.FUNCTIONAL_DISRUPTION, RescueAmenability.LOW
        if rescue_prob >= HIGH_CONFIDENCE_THRESHOLD:
            return MutationClass.MISFOLDING, RescueAmenability.HIGH
        elif rescue_prob >= RESCUE_PROBABILITY_THRESHOLD:
            return MutationClass.MISFOLDING, RescueAmenability.MODERATE
        else:
            return MutationClass.UNKNOWN, RescueAmenability.LOW

    def _classify_with_model(
        self,
        structure: ProteinStructure,
        mutation: Mutation,
        ddg: float,
        esm2_emb: np.ndarray,
    ) -> float:
        """Run the rescue classifier if loaded, else use heuristic."""
        if self.rescue_classifier is None:
            return self._heuristic_rescue_prob(ddg)

        import torch

        try:
            pos_0 = mutation.position - 1
            ca = structure.ca_coords
            mut_ca = ca[pos_0]
            dists = np.linalg.norm(ca - mut_ca, axis=-1)
            neighbor_mask = (dists < CONTACT_CA_CUTOFF) & structure.residue_mask
            indices = np.where(neighbor_mask)[0]

            n_nodes = max(len(indices), 1)
            node_feats = np.zeros((n_nodes, GNN_NODE_DIM), dtype=np.float32)
            for i, ri in enumerate(indices):
                rt = int(structure.residue_types[ri])
                if rt < 20:
                    node_feats[i, rt] = 1.0
                node_feats[i, 20] = structure.bfactors[ri] / 100.0
                node_feats[i, 22] = float(ri == pos_0)

            thermo = np.zeros(THERMO_FEAT_DIM, dtype=np.float32)
            thermo[0] = float(np.clip(ddg / 10.0, -2.0, 2.0))

            with torch.no_grad():
                nf = torch.tensor(node_feats).unsqueeze(0).to(self.device)
                ei = torch.zeros((2, 1), dtype=torch.long).to(self.device)
                ef = torch.zeros((1, GNN_EDGE_DIM)).to(self.device)
                esm2 = torch.tensor(esm2_emb).unsqueeze(0).to(self.device)
                th = torch.tensor(thermo).unsqueeze(0).to(self.device)
                ev = torch.zeros(1, EVO_FEAT_DIM).to(self.device)
                logit = self.rescue_classifier(nf, ei, ef, esm2, th, ev)
                prob = float(torch.sigmoid(logit).item())
            return prob
        except Exception as e:
            logger.debug(f"Classifier forward pass failed: {e}")
            return self._heuristic_rescue_prob(ddg)

    @staticmethod
    def _heuristic_rescue_prob(ddg: float) -> float:
        """Sigmoid heuristic rescue probability from ΔΔG."""
        import math
        if ddg < DDG_DESTABILIZING_THRESHOLD or ddg > DDG_SEVERELY_UNSTABLE:
            return 0.1
        return float(1.0 / (1.0 + math.exp(-(ddg - 2.5))))

    def _get_esm2_embedding(self, sequence: str) -> np.ndarray:
        """ESM-2 150M mean pool representation layer 30."""
        try:
            import torch
            import esm

            if not hasattr(self, "_esm2_model"):
                self._esm2_model, self._esm2_alphabet = esm.pretrained.esm2_t30_150M_UR50D()
                self._esm2_model = self._esm2_model.to(self.device)
                self._esm2_model.eval()

            seq = sequence[:ESM1V_MAX_SEQ_LEN]
            batch_converter = self._esm2_alphabet.get_batch_converter()
            _, _, tokens = batch_converter([("protein", seq)])
            tokens = tokens.to(self.device)

            with torch.no_grad():
                results = self._esm2_model(tokens, repr_layers=[30], return_contacts=False)
                emb = results["representations"][30][0, 1:-1].mean(dim=0).cpu().numpy()
            return emb.astype(np.float32)
        except Exception as e:
            logger.debug(f"ESM-2 embedding failed: {e}")
            return np.zeros(ESM2_EMBED_DIM, dtype=np.float32)

    def _save_result(self, result: REFOLDResult, output_dir: Path) -> None:
        from refold.utils.io import save_result_json, save_molecules_csv
        output_dir.mkdir(parents=True, exist_ok=True)
        mut = result.mutation
        stem = f"{mut.uniprot_id}_{mut.position}_{mut.wildtype_aa}{mut.mutant_aa}"
        save_result_json(result.to_dict(), output_dir / f"{stem}_result.json")
        if result.top_candidates:
            save_molecules_csv(
                [m.to_dict() for m in result.top_candidates],
                output_dir / f"{stem}_candidates.csv",
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="REFOLD single mutation pipeline")
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--position", type=int, required=True)
    parser.add_argument("--wildtype", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--gene", default="")
    parser.add_argument("--disease", default="")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-conformations", type=int, default=N_CONFORMATIONS)
    parser.add_argument("--n-molecules", type=int, default=N_MOLECULES_PER_POCKET)
    parser.add_argument("--device", default=None)
    args = parser.parse_args()

    from refold.utils.logging import setup_logging
    setup_logging()

    device = None
    if args.device:
        import torch
        device = torch.device(args.device)

    mutation = Mutation(
        uniprot_id=args.uniprot,
        position=args.position,
        wildtype_aa=args.wildtype.upper(),
        mutant_aa=args.mutant.upper(),
        gene_name=args.gene,
        disease=args.disease,
    )

    pipeline = REFOLDPipeline(
        device=device,
        n_conformations=args.n_conformations,
        n_molecules_per_pocket=args.n_molecules,
    )
    result = pipeline.run(mutation, output_dir=args.output_dir)

    print(f"\nREFOLD Result for {mutation.hgvs} ({mutation.gene_name or mutation.uniprot_id})")
    print(f"  Mutation class:      {result.mutation_class.value}")
    print(f"  Rescue amenability:  {result.rescue_amenability.value}")
    print(f"  Rescue probability:  {result.rescue_probability:.3f}")
    print(f"  ΔΔG predicted:       {result.ddg_predicted:.2f} kcal/mol")
    print(f"  Pockets detected:    {result.n_pockets_detected}")
    print(f"  Molecules generated: {result.n_molecules_generated}")
    if result.top_candidates:
        best = result.top_candidates[0]
        print(f"\n  Top candidate:")
        print(f"    SMILES: {best.smiles}")
        print(f"    QED:    {best.qed_score:.3f}")
        print(f"    SA:     {best.sa_score:.2f}")
        print(f"    RP:     {best.rescue_probability:.3f}")
    print(f"\n  Runtime: {result.runtime_seconds:.1f}s")


if __name__ == "__main__":
    main()
