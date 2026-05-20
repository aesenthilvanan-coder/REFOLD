"""Solvent-accessible surface area computation."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.constants import MAX_ASA, SURFACE_PROBE_RADIUS, STANDARD_AAS
from refold.types import ProteinStructure

logger = logging.getLogger(__name__)

# Fibonacci sphere for SASA estimation
_N_SPHERE_POINTS = 960


def _fibonacci_sphere(n: int) -> np.ndarray:
    """Generate *n* approximately evenly-spaced unit vectors on a sphere."""
    golden = (1 + np.sqrt(5)) / 2
    i = np.arange(n)
    theta = np.arccos(1 - 2 * (i + 0.5) / n)
    phi = 2 * np.pi * i / golden
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.stack([x, y, z], axis=1).astype(np.float32)


_SPHERE_POINTS = _fibonacci_sphere(_N_SPHERE_POINTS)

# Van der Waals radii (Å) for common heavy atoms
_VDW_RADII: dict[str, float] = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "Cl": 1.75, "Br": 1.85,
    "I": 1.98, "H": 1.20,
}


def compute_rel_asa(structure: ProteinStructure) -> np.ndarray:
    """
    Estimate per-residue relative solvent-accessible surface area (0-1)
    using a simplified shrake-rupley approach on Cα only.

    This is a fast approximation. For publication-quality SASA use DSSP or FreeSASA.
    """
    n = structure.n_residues
    ca = structure.ca_coords.copy()  # [N, 3]
    probe = SURFACE_PROBE_RADIUS
    r_ca = _VDW_RADII["C"] + probe  # effective radius per Cα

    rel_asa = np.zeros(n, dtype=np.float32)

    for i in range(n):
        if np.any(np.isnan(ca[i])):
            rel_asa[i] = np.nan
            continue

        # Sphere points around residue i
        test_pts = ca[i] + r_ca * _SPHERE_POINTS  # [K, 3]

        # Find neighbors within 2*r_ca + probe
        cutoff = 2 * r_ca + probe
        diff = ca - ca[i]  # [N, 3]
        dist_sq = np.sum(diff ** 2, axis=-1)
        neighbors = np.where((dist_sq < cutoff ** 2) & (dist_sq > 0))[0]

        accessible = np.ones(len(_SPHERE_POINTS), dtype=bool)
        for j in neighbors:
            if np.any(np.isnan(ca[j])):
                continue
            d = np.sqrt(np.sum((test_pts - ca[j]) ** 2, axis=-1))
            accessible &= d > r_ca

        asa_i = 4 * np.pi * r_ca ** 2 * accessible.sum() / len(_SPHERE_POINTS)
        aa = structure.sequence[i] if i < len(structure.sequence) else "G"
        max_asa = MAX_ASA.get(aa, MAX_ASA["G"])
        rel_asa[i] = float(np.clip(asa_i / max_asa, 0.0, 1.0))

    return rel_asa


def get_surface_residues(
    structure: ProteinStructure,
    rel_asa_threshold: float = 0.25,
    rel_asa: Optional[np.ndarray] = None,
) -> list[int]:
    """Return 0-based indices of solvent-exposed residues."""
    if rel_asa is None:
        rel_asa = compute_rel_asa(structure)
    return [i for i, v in enumerate(rel_asa) if not np.isnan(v) and v >= rel_asa_threshold]
