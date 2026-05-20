"""Rescue probability computation and molecule ranking."""

from __future__ import annotations

import math
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from refold.types import GeneratedMolecule, Pocket

logger = logging.getLogger(__name__)


def compute_final_rescue_probability(
    molecule: "GeneratedMolecule",
    pocket: "Pocket",
    rescue_amenability_prob: float,
    ddg_protein: float,
) -> float:
    """Compute final rescue probability as weighted combination of signals.

    Component weights:
        s_rescue        0.30  — classifier rescue amenability
        s_affinity      0.25  — binding affinity (sigmoid at -5 kcal/mol)
        s_druggability  0.15  — pocket druggability score
        s_detection     0.10  — pocket detection frequency
        s_transient     0.10  — transient pocket bonus
        s_druglike      0.05  — drug-likeness (QED)
        s_severity      0.05  — mutation severity (ΔΔG contribution)

    Returns: probability in [0, 1].
    """
    s_rescue = float(rescue_amenability_prob)

    affinity = molecule.predicted_affinity_kcal or 0.0
    s_affinity = 1.0 / (1.0 + math.exp(affinity + 5.0))

    s_druggability = float(pocket.druggability_score)

    s_detection = float(pocket.detection_frequency)

    s_transient = 1.0 if pocket.is_transient else 0.5

    s_druglike = float(molecule.qed_score) if molecule.qed_score else 0.0

    s_severity = min(abs(ddg_protein) / 5.0, 1.0) if ddg_protein > 0 else 0.0

    rescue_prob = (
        0.30 * s_rescue
        + 0.25 * s_affinity
        + 0.15 * s_druggability
        + 0.10 * s_detection
        + 0.10 * s_transient
        + 0.05 * s_druglike
        + 0.05 * s_severity
    )

    return float(max(0.0, min(1.0, rescue_prob)))


def rank_molecules_by_rescue_probability(
    molecules: list["GeneratedMolecule"],
) -> list["GeneratedMolecule"]:
    """Sort molecules by rescue probability descending."""
    ranked = sorted(molecules, key=lambda m: m.rescue_probability, reverse=True)
    for i, mol in enumerate(ranked):
        mol.rank = i + 1
    return ranked
