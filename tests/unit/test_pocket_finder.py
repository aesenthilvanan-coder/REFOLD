"""Unit tests for pocket finder."""

import numpy as np
import pytest

from refold.structure.pocket_finder import (
    _cluster_pockets_across_conformations,
    _compute_hydrophobicity,
)


def test_cluster_empty_pocket_sets():
    clusters = _cluster_pockets_across_conformations([[], [], []])
    assert clusters == []


def test_cluster_single_pocket():
    center = np.array([10.0, 5.0, 3.0], dtype=np.float32)
    spheres = np.tile(center, (5, 1))
    pocket = {"alpha_spheres": spheres, "volume": 500.0, "druggability": 0.7}
    clusters = _cluster_pockets_across_conformations([[pocket], [pocket], []])
    assert len(clusters) == 1
    assert clusters[0]["detection_frequency"] == pytest.approx(2 / 3)


def test_cluster_merges_nearby():
    center1 = np.array([10.0, 5.0, 3.0], dtype=np.float32)
    center2 = np.array([11.0, 5.5, 3.0], dtype=np.float32)
    s1 = np.tile(center1, (5, 1))
    s2 = np.tile(center2, (5, 1))
    p1 = {"alpha_spheres": s1, "volume": 500.0, "druggability": 0.7}
    p2 = {"alpha_spheres": s2, "volume": 510.0, "druggability": 0.72}
    clusters = _cluster_pockets_across_conformations([[p1], [p2]])
    assert len(clusters) == 1


def test_compute_hydrophobicity(dummy_structure):
    score = _compute_hydrophobicity(dummy_structure, list(range(10)))
    assert 0.0 <= score <= 1.0
