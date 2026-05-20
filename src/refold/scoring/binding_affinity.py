"""Predicted binding affinity scoring for pocket-molecule pairs."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.types import Pocket, GeneratedMolecule

logger = logging.getLogger(__name__)


def estimate_binding_affinity_heuristic(
    molecule: GeneratedMolecule,
    pocket: Pocket,
) -> float:
    """
    Fast heuristic binding affinity estimate (kcal/mol, negative = favorable).

    Uses a simple scoring function based on:
    - QED (proxy for drug-like binding)
    - Pocket druggability
    - Hydrophobic matching
    - MW penalty

    Returns predicted ΔG in kcal/mol.
    """
    # Base affinity from pocket druggability
    base = -2.0 * pocket.druggability_score

    # Hydrophobic complement
    from refold.constants import MAX_LOGP
    logp_norm = float(np.clip(molecule.logp / MAX_LOGP, 0.0, 1.0))
    hydrophobic_bonus = -1.5 * pocket.hydrophobicity * logp_norm

    # QED quality bonus
    qed_bonus = -2.0 * molecule.qed_score

    # Transient pocket bonus (harder to compete → better selectivity)
    from refold.types import PocketType
    transient_bonus = -0.5 if pocket.pocket_type == PocketType.TRANSIENT_MISFOLDING else 0.0

    # MW entropy penalty (too large → entropic cost)
    from refold.constants import MAX_MW
    mw_penalty = 0.5 * float(np.clip((molecule.mw - 300) / (MAX_MW - 300), 0.0, 1.0))

    affinity = base + hydrophobic_bonus + qed_bonus + transient_bonus + mw_penalty
    return float(np.clip(affinity, -15.0, 0.0))


def score_molecule_pocket_affinity(
    smiles: str,
    pocket: Pocket,
    method: str = "heuristic",
) -> float:
    """
    Score binding affinity for SMILES + pocket pair.

    Args:
        smiles: SMILES string
        pocket: Pocket object
        method: "heuristic" (default) or future ML methods

    Returns:
        Predicted ΔG in kcal/mol.
    """
    if method == "heuristic":
        # Build a minimal molecule object for scoring
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return -3.0
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            qed = 0.5
            try:
                from rdkit.Chem import QED
                qed = QED.qed(mol)
            except Exception:
                pass

            dummy = GeneratedMolecule(smiles=smiles, pocket_id=pocket.pocket_id, mw=mw, logp=logp, qed_score=qed)
            return estimate_binding_affinity_heuristic(dummy, pocket)
        except Exception as e:
            logger.debug(f"Affinity scoring failed: {e}")
            return -3.0
    else:
        raise ValueError(f"Unknown affinity scoring method: {method}")


def affinity_to_ki(affinity_kcal: float, temperature_k: float = 298.15) -> float:
    """Convert predicted ΔG (kcal/mol) to Ki (μM)."""
    from refold.utils.chemistry import kcal_to_um
    return kcal_to_um(affinity_kcal, temperature_k)
