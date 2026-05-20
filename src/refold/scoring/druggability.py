"""Pocket druggability scoring."""

from __future__ import annotations

import logging
import numpy as np

from refold.types import Pocket, PocketType
from refold.constants import DRUGGABILITY_THRESHOLD, MIN_POCKET_VOLUME, MAX_POCKET_VOLUME

logger = logging.getLogger(__name__)


def compute_druggability_score(
    volume: float,
    fpocket_score: float,
    hydrophobicity: float,
    detection_freq: float,
    pocket_type: PocketType,
    n_residues: int = 0,
) -> float:
    """
    Compute composite druggability score [0, 1].

    Volume component peaks in the 400-900 Å³ sweet spot.
    Transient misfolding pockets receive a bonus.
    """
    # Volume score — favor 400-900 Å³ sweet spot
    if volume < MIN_POCKET_VOLUME:
        vol_score = 0.0
    elif volume > MAX_POCKET_VOLUME:
        vol_score = max(0.0, 1.0 - (volume - MAX_POCKET_VOLUME) / MAX_POCKET_VOLUME)
    else:
        center = (MIN_POCKET_VOLUME + MAX_POCKET_VOLUME) / 2.0
        width = (MAX_POCKET_VOLUME - MIN_POCKET_VOLUME) / 2.0
        vol_score = 1.0 - abs(volume - center) / (width + 1.0)

    # Residue count score
    res_score = float(np.clip(n_residues / 20.0, 0.0, 1.0)) if n_residues > 0 else 0.5

    # Detection reliability
    freq_score = float(np.clip(detection_freq, 0.0, 1.0))

    transient_bonus = 0.1 if pocket_type == PocketType.TRANSIENT_MISFOLDING else 0.0

    score = (
        0.35 * float(np.clip(fpocket_score, 0.0, 1.0))
        + 0.25 * float(np.clip(hydrophobicity, 0.0, 1.0))
        + 0.20 * float(np.clip(vol_score, 0.0, 1.0))
        + 0.10 * freq_score
        + 0.10 * res_score
        + transient_bonus
    )
    return float(np.clip(score, 0.0, 1.0))


def classify_druggability(score: float) -> str:
    """Map druggability score to human-readable tier."""
    if score >= 0.75:
        return "highly_druggable"
    elif score >= DRUGGABILITY_THRESHOLD:
        return "druggable"
    elif score >= 0.30:
        return "moderately_druggable"
    else:
        return "undruggable"


def score_all_pockets(pockets: list[Pocket]) -> list[tuple[Pocket, float, str]]:
    """Return [(pocket, score, tier)] sorted by score descending."""
    results = []
    for p in pockets:
        score = compute_druggability_score(
            volume=p.volume,
            fpocket_score=p.druggability_score,
            hydrophobicity=p.hydrophobicity,
            detection_freq=p.detection_frequency,
            pocket_type=p.pocket_type,
            n_residues=len(p.residue_indices),
        )
        tier = classify_druggability(score)
        results.append((p, score, tier))
    results.sort(key=lambda x: x[1], reverse=True)
    return results
