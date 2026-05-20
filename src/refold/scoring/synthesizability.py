"""Synthesizability scoring for generated molecules."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.constants import MAX_SA_SCORE

logger = logging.getLogger(__name__)


def compute_sa_score(smiles: str) -> Optional[float]:
    """Compute RDKit SA score (1=easy, 10=hard). Returns None on failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import RDConfig
        import sys
        import os
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return float(sascorer.calculateScore(mol))
    except Exception:
        return _fallback_sa_score(smiles)


def _fallback_sa_score(smiles: str) -> Optional[float]:
    """Heuristic SA score fallback based on ring count and complexity."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Penalty for ring complexity and stereocenters
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        score = 2.0 + 0.5 * n_rings + 0.7 * n_stereo + 1.0 * n_spiro
        return float(np.clip(score, 1.0, 10.0))
    except Exception as e:
        logger.debug(f"Fallback SA score failed: {e}")
        return None


def compute_synthetic_accessibility(smiles: str) -> float:
    """
    Return synthesizability score [0, 1] where 1 = easy to synthesize.

    Uses SA score normalized and inverted.
    """
    sa = compute_sa_score(smiles)
    if sa is None:
        return 0.5
    # SA score: 1 (easy) → 10 (hard). Normalize to [0, 1] inverted.
    return float(np.clip(1.0 - (sa - 1.0) / 9.0, 0.0, 1.0))


def passes_synthesizability_filter(smiles: str, max_sa_score: float = MAX_SA_SCORE) -> bool:
    """Return True if SA score ≤ threshold."""
    sa = compute_sa_score(smiles)
    if sa is None:
        return True  # benefit of the doubt
    return sa <= max_sa_score


def estimate_retrosynthetic_steps(smiles: str) -> int:
    """Rough estimate of retrosynthetic steps based on molecular complexity."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 5
        n_heavy = mol.GetNumHeavyAtoms()
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        steps = 1 + n_rings + max(0, (n_heavy - 10) // 5)
        return int(np.clip(steps, 1, 15))
    except Exception:
        return 5
