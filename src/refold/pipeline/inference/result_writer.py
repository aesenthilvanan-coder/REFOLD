"""Write REFOLD results to various output formats."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from refold.types import REFOLDResult

logger = logging.getLogger(__name__)


def write_json(result: REFOLDResult, path: Path) -> None:
    """Write single result to JSON."""
    from refold.utils.io import NumpyEncoder
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(result.to_dict(), f, cls=NumpyEncoder, indent=2)


def write_results_to_parquet(results: list[REFOLDResult], path: Path) -> None:
    """Flatten results to tabular DataFrame and write parquet."""
    rows = []
    for r in results:
        row = {
            "uniprot_id": r.mutation.uniprot_id,
            "gene_name": r.mutation.gene_name,
            "position": r.mutation.position,
            "wildtype_aa": r.mutation.wildtype_aa,
            "mutant_aa": r.mutation.mutant_aa,
            "hgvs": r.mutation.hgvs,
            "disease": r.mutation.disease,
            "clinvar_id": r.mutation.clinvar_id,
            "mutation_class": r.mutation_class.value,
            "rescue_amenability": r.rescue_amenability.value,
            "rescue_amenability_prob": r.rescue_amenability_prob,
            "rescue_probability": r.rescue_probability,
            "ddg_predicted": r.ddg_predicted,
            "n_pockets_detected": r.n_pockets_detected,
            "n_molecules_generated": r.n_molecules_generated,
            "n_molecules_passing_filters": r.n_molecules_passing_filters,
            "runtime_seconds": r.runtime_seconds,
            "error_message": r.error_message,
        }
        # Top candidate info
        if r.top_candidates:
            top = r.top_candidates[0]
            row["top_smiles"] = top.smiles
            row["top_qed"] = top.qed_score
            row["top_sa_score"] = top.sa_score
            row["top_rescue_prob"] = top.rescue_probability
            row["top_affinity_kcal"] = top.predicted_affinity_kcal
        else:
            row["top_smiles"] = None
            row["top_qed"] = None
            row["top_sa_score"] = None
            row["top_rescue_prob"] = None
            row["top_affinity_kcal"] = None

        rows.append(row)

    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    logger.info(f"Wrote {len(df)} results to {path}")


def write_csv_summary(results: list[REFOLDResult], path: Path) -> None:
    """Write human-readable CSV summary."""
    rows = []
    for r in results:
        top_smiles = r.top_candidates[0].smiles if r.top_candidates else ""
        top_qed = f"{r.top_candidates[0].qed_score:.3f}" if r.top_candidates else ""
        rows.append({
            "mutation": str(r.mutation),
            "rescue_amenability": r.rescue_amenability.value,
            "rescue_prob": f"{r.rescue_probability:.3f}",
            "ddg": f"{r.ddg_predicted:.2f}",
            "n_pockets": r.n_pockets_detected,
            "n_candidates": r.n_molecules_passing_filters,
            "top_smiles": top_smiles,
            "top_qed": top_qed,
            "runtime_s": f"{r.runtime_seconds:.1f}",
        })
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    logger.info(f"CSV summary written to {path}")


def write_rescue_candidates_sdf(results: list[REFOLDResult], path: Path) -> None:
    """Write top rescue candidate molecules to SDF format."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, SDWriter
    except ImportError:
        logger.warning("RDKit required for SDF output")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    writer = SDWriter(str(path))
    n_written = 0

    for r in results:
        for mol_obj in r.top_candidates:
            try:
                mol = Chem.MolFromSmiles(mol_obj.smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)
                mol = Chem.RemoveHs(mol)
                mol.SetProp("_Name", f"{r.mutation.to_key()}_rank{mol_obj.rank}")
                mol.SetProp("SMILES", mol_obj.smiles)
                mol.SetProp("rescue_probability", f"{mol_obj.rescue_probability:.4f}")
                mol.SetProp("qed_score", f"{mol_obj.qed_score:.4f}")
                mol.SetProp("sa_score", f"{mol_obj.sa_score:.4f}")
                mol.SetProp("pocket_id", mol_obj.pocket_id)
                mol.SetProp("uniprot_id", r.mutation.uniprot_id)
                writer.write(mol)
                n_written += 1
            except Exception as e:
                logger.debug(f"SDF write failed for {mol_obj.smiles}: {e}")

    writer.close()
    logger.info(f"Wrote {n_written} molecules to {path}")
