"""Pocket encoder for conditioning molecule generation."""

from __future__ import annotations

import numpy as np


def build_pocket_node_features(pocket) -> np.ndarray:
    """
    Build [n_spheres, POCKET_NODE_FEAT_DIM] node feature matrix for alpha spheres.

    Features per sphere:
      [0:3]  normalized XYZ coordinates (divided by 50)
      [3]    sphere radius placeholder (0.5)
      [4]    druggability score
      [5]    hydrophobicity
      [6]    detection frequency
      [7]    is_transient flag
      [8:25] reserved zeros
    """
    from refold.constants import POCKET_NODE_FEAT_DIM
    from refold.types import PocketType

    spheres = pocket.alpha_sphere_coords  # [K, 3]
    k = len(spheres)

    feats = np.zeros((k, POCKET_NODE_FEAT_DIM), dtype=np.float32)
    center = pocket.center

    for i, coord in enumerate(spheres):
        feats[i, 0:3] = (coord - center) / 50.0
        feats[i, 3] = 0.5
        feats[i, 4] = pocket.druggability_score
        feats[i, 5] = pocket.hydrophobicity
        feats[i, 6] = pocket.detection_frequency
        feats[i, 7] = float(pocket.pocket_type == PocketType.TRANSIENT_MISFOLDING)

    return feats


def build_pocket_scalar_features(pocket) -> np.ndarray:
    """
    Build [1, POCKET_SCALAR_DIM] global pocket scalar features.

    Features:
      [0] volume / 1000
      [1] druggability score
      [2] hydrophobicity / 5
      [3] detection frequency
      [4:7] center coordinates / 50
      [7] is_transient
    """
    from refold.constants import POCKET_SCALAR_DIM
    from refold.types import PocketType

    center = pocket.center
    feats = np.array([
        pocket.volume / 1000.0,
        pocket.druggability_score,
        pocket.hydrophobicity / 5.0,
        pocket.detection_frequency,
        center[0] / 50.0,
        center[1] / 50.0,
        center[2] / 50.0,
        float(pocket.pocket_type == PocketType.TRANSIENT_MISFOLDING),
    ], dtype=np.float32)

    assert len(feats) == POCKET_SCALAR_DIM
    return feats.reshape(1, -1)
