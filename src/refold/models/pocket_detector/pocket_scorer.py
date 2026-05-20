"""Score and rank detected pockets for druggability."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.constants import (
    DRUGGABILITY_THRESHOLD, MIN_POCKET_VOLUME, MAX_POCKET_VOLUME,
    TRANSIENT_POCKET_FREQ_THRESHOLD,
)
from refold.types import Pocket, PocketType

logger = logging.getLogger(__name__)


def score_pocket_druggability(
    volume: float,
    druggability_score: float,
    hydrophobicity: float,
    detection_freq: float,
    pocket_type: PocketType,
) -> float:
    """
    Compute a composite druggability score [0, 1].

    Weighted combination of:
    - fpocket druggability score (0.4)
    - hydrophobicity (0.25)
    - volume score — Goldilocks range (0.20)
    - detection frequency (0.10)
    - bonus for transient pocket (0.05 extra weight)
    """
    # Volume score: peaks in 300-2000 Å³ range
    vol_score = 0.0
    if MIN_POCKET_VOLUME <= volume <= MAX_POCKET_VOLUME:
        vol_norm = (volume - MIN_POCKET_VOLUME) / (MAX_POCKET_VOLUME - MIN_POCKET_VOLUME)
        vol_score = 1.0 - abs(vol_norm - 0.4) / 0.6  # peak near 600-800 Å³

    transient_bonus = 0.1 if pocket_type == PocketType.TRANSIENT_MISFOLDING else 0.0

    score = (
        0.40 * float(np.clip(druggability_score, 0.0, 1.0))
        + 0.25 * float(np.clip(hydrophobicity, 0.0, 1.0))
        + 0.20 * float(np.clip(vol_score, 0.0, 1.0))
        + 0.10 * float(np.clip(detection_freq, 0.0, 1.0))
        + transient_bonus
    )
    return float(np.clip(score, 0.0, 1.0))


def rank_pockets(pockets: list[Pocket]) -> list[Pocket]:
    """Return pockets sorted by composite druggability score descending."""
    def _sort_key(p: Pocket) -> float:
        composite = score_pocket_druggability(
            p.volume, p.druggability_score, p.hydrophobicity,
            p.detection_frequency, p.pocket_type,
        )
        # Transient misfolding pockets are always preferred
        return composite + (0.2 if p.pocket_type == PocketType.TRANSIENT_MISFOLDING else 0.0)

    return sorted(pockets, key=_sort_key, reverse=True)


def filter_pockets(
    pockets: list[Pocket],
    min_volume: float = MIN_POCKET_VOLUME,
    max_volume: float = MAX_POCKET_VOLUME,
    min_druggability: float = DRUGGABILITY_THRESHOLD,
    require_transient: bool = False,
) -> list[Pocket]:
    """Filter pockets by volume and druggability thresholds."""
    result = []
    for p in pockets:
        if p.volume < min_volume or p.volume > max_volume:
            continue
        if p.druggability_score < min_druggability:
            continue
        if require_transient and p.pocket_type != PocketType.TRANSIENT_MISFOLDING:
            continue
        result.append(p)
    return result


def select_top_pockets(
    pockets: list[Pocket],
    n_top: int = 3,
    prefer_transient: bool = True,
) -> list[Pocket]:
    """Select the top-n pockets, optionally prioritising transient pockets."""
    if not pockets:
        return []

    ranked = rank_pockets(pockets)

    if prefer_transient:
        transient = [p for p in ranked if p.pocket_type == PocketType.TRANSIENT_MISFOLDING]
        others = [p for p in ranked if p.pocket_type != PocketType.TRANSIENT_MISFOLDING]
        combined = transient + others
    else:
        combined = ranked

    return combined[:n_top]
