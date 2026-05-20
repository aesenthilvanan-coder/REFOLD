"""Mutation application and ΔΔG prediction.

ESM-1v: masked marginal probability proxy.
GNN: structure-based ΔΔG.
Ensemble: weighted combination (0.4 ESM-1v + 0.6 GNN).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from refold.constants import (
    CA_IDX, N_IDX, C_IDX, O_IDX, CB_IDX,
    N_ATOM_TYPES, AA_ONE_TO_THREE, AA_THREE_TO_ONE,
    STANDARD_AAS, ESM1V_MAX_SEQ_LEN,
    CONTACT_CA_CUTOFF, GNN_NODE_DIM, GNN_EDGE_DIM,
)

if TYPE_CHECKING:
    from refold.types import ProteinStructure, Mutation, MutantStructure

logger = logging.getLogger(__name__)


class ESM1vStabilityPredictor:
    """ΔΔG proxy via ESM-1v masked marginal probability.

    ΔΔG = (log P(wt|context) - log P(mut|context)) × 2.5
    Positive = destabilizing.
    """

    def __init__(self, device: "torch.device | None" = None):
        self._model = None
        self._alphabet = None
        self._device = device

    def _load_model(self) -> None:
        if self._model is not None:
            return
        import torch
        import esm

        if self._device is None:
            from refold.utils.device import get_device
            self._device = get_device()

        self._model, self._alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
        self._model = self._model.to(self._device)
        self._model.eval()
        logger.info("ESM-1v model loaded")

    def predict_ddg(
        self,
        sequence: str,
        position: int,
        wildtype_aa: str,
        mutant_aa: str,
    ) -> float:
        """Predict ΔΔG for a single mutation using masked marginal probability.

        Args:
            sequence: Full protein sequence.
            position: 1-based residue position.
            wildtype_aa: Single-letter wildtype AA.
            mutant_aa: Single-letter mutant AA.

        Returns:
            ΔΔG in kcal/mol (positive = destabilizing).
        """
        import torch

        try:
            self._load_model()
        except Exception as e:
            logger.warning(f"ESM-1v not available: {e}")
            return 0.0

        seq = sequence[:ESM1V_MAX_SEQ_LEN]
        pos_0 = position - 1
        if pos_0 < 0 or pos_0 >= len(seq):
            return 0.0

        batch_converter = self._alphabet.get_batch_converter()
        masked_seq = seq[:pos_0] + self._alphabet.mask_token + seq[pos_0 + 1:]
        data = [("protein", masked_seq)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(self._device)

        wt_idx = self._alphabet.get_idx(wildtype_aa)
        mt_idx = self._alphabet.get_idx(mutant_aa)
        mask_pos = pos_0 + 1  # account for BOS token

        with torch.no_grad():
            logits = self._model(tokens)["logits"]
            log_probs = torch.log_softmax(logits[0, mask_pos], dim=-1)
            log_p_wt = log_probs[wt_idx].item()
            log_p_mt = log_probs[mt_idx].item()

        ddg = (log_p_wt - log_p_mt) * 2.5
        return float(ddg)


class StructureBasedDDGPredictor:
    """GNN-based ΔΔG prediction from local structural neighborhood."""

    def __init__(self, checkpoint_path: str | None = None, device: "torch.device | None" = None):
        self._checkpoint = checkpoint_path
        self._model = None
        self._device = device

    def _load_model(self) -> None:
        if self._model is not None:
            return
        import torch
        from refold.constants import CHECKPOINT_DIR

        if self._device is None:
            from refold.utils.device import get_device
            self._device = get_device()

        ckpt_path = self._checkpoint or str(CHECKPOINT_DIR / "rescue_classifier" / "best.pt")

        try:
            from refold.models.rescue_classifier.model import RescueClassifier
            ckpt = torch.load(ckpt_path, map_location=self._device)
            self._model = RescueClassifier()
            self._model.load_state_dict(ckpt.get("model_state_dict", ckpt), strict=False)
            self._model = self._model.to(self._device)
            self._model.eval()
        except Exception as e:
            logger.debug(f"Structure-based predictor checkpoint not loaded: {e}")
            self._model = None

    def predict_ddg(
        self,
        structure: "ProteinStructure",
        position: int,
        wildtype_aa: str,
        mutant_aa: str,
    ) -> float:
        """Predict ΔΔG using GNN on local structural neighborhood."""
        self._load_model()

        if self._model is None:
            return 0.0

        import torch

        pos_0 = position - 1
        ca = structure.ca_coords
        mut_ca = ca[pos_0]
        dists = np.linalg.norm(ca - mut_ca, axis=-1)
        neighbor_mask = dists < CONTACT_CA_CUTOFF

        neighbor_indices = np.where(neighbor_mask)[0]
        if len(neighbor_indices) == 0:
            return 0.0

        node_feats = np.zeros((len(neighbor_indices), GNN_NODE_DIM), dtype=np.float32)
        for i, ri in enumerate(neighbor_indices):
            rt = int(structure.residue_types[ri])
            if rt < 20:
                node_feats[i, rt] = 1.0
            node_feats[i, 20] = structure.bfactors[ri] / 100.0
            node_feats[i, 22] = float(ri == pos_0)

        src, dst, edge_feats = [], [], []
        for i, ri in enumerate(neighbor_indices):
            for j, rj in enumerate(neighbor_indices):
                if i == j:
                    continue
                d = np.linalg.norm(ca[ri] - ca[rj])
                if d < CONTACT_CA_CUTOFF:
                    src.append(i)
                    dst.append(j)
                    ef = np.zeros(GNN_EDGE_DIM, dtype=np.float32)
                    ef[0] = d / CONTACT_CA_CUTOFF
                    ef[1] = abs(ri - rj) / len(ca)
                    edge_feats.append(ef)

        if not src:
            return 0.0

        try:
            with torch.no_grad():
                nf = torch.tensor(node_feats).unsqueeze(0).to(self._device)
                ei = torch.tensor([src, dst], dtype=torch.long).to(self._device)
                ef_t = torch.tensor(np.array(edge_feats)).to(self._device)
                esm2 = torch.zeros(1, 480).to(self._device)
                thermo = torch.zeros(1, 16).to(self._device)
                evo = torch.zeros(1, 32).to(self._device)
                out = self._model(nf, ei, ef_t, esm2, thermo, evo)
                return float(out.squeeze().item())
        except Exception as e:
            logger.debug(f"GNN forward pass failed: {e}")
            return 0.0


def predict_ddg_ensemble(
    sequence: str,
    structure: "ProteinStructure",
    mutation: "Mutation",
    esm1v_predictor: ESM1vStabilityPredictor,
    structure_predictor: StructureBasedDDGPredictor,
    esm1v_weight: float = 0.4,
    structure_weight: float = 0.6,
) -> tuple[float, float, float]:
    """Predict ΔΔG using weighted ensemble of ESM-1v and structure-based GNN.

    Returns (ddg_ensemble, ddg_esm1v, ddg_gnn).
    """
    ddg_esm1v = esm1v_predictor.predict_ddg(
        sequence, mutation.position, mutation.wildtype_aa, mutation.mutant_aa
    )
    ddg_gnn = structure_predictor.predict_ddg(
        structure, mutation.position, mutation.wildtype_aa, mutation.mutant_aa
    )

    ddg_ensemble = esm1v_weight * ddg_esm1v + structure_weight * ddg_gnn
    return ddg_ensemble, ddg_esm1v, ddg_gnn


def apply_mutation(
    structure: "ProteinStructure",
    mutation: "Mutation",
) -> "MutantStructure":
    """Apply a point mutation to a ProteinStructure.

    Erases side-chain atoms (sets to NaN), updates residue type index.
    Backbone atoms (N=0, CA=1, C=2, O=4) are preserved.
    """
    from refold.types import MutantStructure
    from refold.constants import AA_TO_IDX

    pos_0 = mutation.position - 1
    n_res = len(structure.sequence)

    if pos_0 < 0 or pos_0 >= n_res:
        raise ValueError(
            f"Mutation position {mutation.position} out of range [1, {n_res}]"
        )

    mutant_coords = structure.coords.copy()
    mutant_residue_types = structure.residue_types.copy()

    # Erase side-chain atoms (keep backbone: N=0, CA=1, C=2, O=4)
    BACKBONE_INDICES = {N_IDX, CA_IDX, C_IDX, O_IDX}
    for atom_idx in range(N_ATOM_TYPES):
        if atom_idx not in BACKBONE_INDICES:
            mutant_coords[pos_0, atom_idx, :] = np.nan

    # Update residue type: AA_TO_IDX maps single-letter → 0-19, unknown → 20
    mt_type_idx = AA_TO_IDX.get(mutation.mutant_aa, 20)
    mutant_residue_types[pos_0] = mt_type_idx

    return MutantStructure(
        wildtype=structure,
        mutation=mutation,
        mutant_coords=mutant_coords,
        mutant_residue_types=mutant_residue_types,
        ddg_predicted=0.0,
        relaxed=False,
    )
