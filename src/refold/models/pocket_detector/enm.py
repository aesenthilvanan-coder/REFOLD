"""Elastic Network Model (ANM) for conformational sampling — exposed as a standalone module."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from scipy import linalg

from refold.constants import (
    ENM_CUTOFF, N_ENM_MODES, ENM_AMPLITUDE_SCALE, ENM_N_SKIP_TRIVIAL,
    N_CONFORMATIONS, CA_IDX,
)
from refold.types import ProteinStructure

logger = logging.getLogger(__name__)


def build_anm_hessian(ca_coords: np.ndarray, cutoff: float = ENM_CUTOFF) -> np.ndarray:
    """
    Construct the 3N×3N ANM Hessian matrix.

    Args:
        ca_coords: [N, 3] Cα coordinates (no NaN)
        cutoff: interaction cutoff in Å

    Returns:
        H: [3N, 3N] symmetric float64 Hessian
    """
    n = len(ca_coords)
    H = np.zeros((3 * n, 3 * n), dtype=np.float64)

    for i in range(n):
        for j in range(i + 1, n):
            diff = ca_coords[i] - ca_coords[j]
            dist2 = float(np.dot(diff, diff))
            if dist2 > cutoff ** 2:
                continue
            # Kirchhoff super-element
            k = np.outer(diff, diff) / dist2  # [3, 3]
            H[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= k
            H[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= k
            H[3 * i:3 * i + 3, 3 * i:3 * i + 3] += k
            H[3 * j:3 * j + 3, 3 * j:3 * j + 3] += k

    return H


def compute_anm_modes(
    ca_coords: np.ndarray,
    n_modes: int = N_ENM_MODES,
    cutoff: float = ENM_CUTOFF,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute ANM normal modes.

    Returns:
        eigenvalues: [n_modes] float64 (non-trivial, ascending)
        eigenvectors: [3N, n_modes] float64
    """
    H = build_anm_hessian(ca_coords, cutoff)
    # Use eigh (symmetric) — faster and guaranteed real
    vals, vecs = linalg.eigh(H)

    # Skip trivial (near-zero) modes
    skip = ENM_N_SKIP_TRIVIAL
    vals = vals[skip:skip + n_modes]
    vecs = vecs[:, skip:skip + n_modes]

    return vals, vecs


def sample_conformations(
    structure: ProteinStructure,
    n_conformations: int = N_CONFORMATIONS,
    n_modes: int = N_ENM_MODES,
    amplitude_scale: float = ENM_AMPLITUDE_SCALE,
    seed: Optional[int] = None,
) -> list[np.ndarray]:
    """
    Generate an ensemble of conformations by ANM displacement.

    Returns list of n_conformations coord arrays [N_res, 37, 3], first being the original.
    """
    rng = np.random.default_rng(seed)

    ca = structure.ca_coords.copy()  # [N, 3]
    # Use only resolved residues
    valid = ~np.any(np.isnan(ca), axis=-1)
    valid_indices = np.where(valid)[0]

    conformations = [structure.coords.copy()]

    if len(valid_indices) < 4 or n_conformations <= 1:
        return conformations + [structure.coords.copy() for _ in range(n_conformations - 1)]

    ca_valid = ca[valid_indices]
    try:
        eigenvalues, eigenvectors = compute_anm_modes(ca_valid, n_modes=n_modes)
    except Exception as e:
        logger.warning(f"ANM mode computation failed: {e}; returning unperturbed conformations")
        return conformations + [structure.coords.copy() for _ in range(n_conformations - 1)]

    # pLDDT-based amplitude scaling: lower confidence → larger displacement
    plddt_valid = structure.bfactors[valid_indices] / 100.0  # [K]
    amplitude_factors = amplitude_scale * (2.0 - np.clip(plddt_valid, 0.3, 1.0))  # [K]

    n_valid = len(valid_indices)

    for _ in range(n_conformations - 1):
        # Random linear combination of modes, amplitude ∝ 1/√eigenvalue
        mode_weights = rng.standard_normal(len(eigenvalues))
        mode_scales = amplitude_scale / (np.sqrt(eigenvalues) + 1e-8)
        displacement_flat = eigenvectors @ (mode_weights * mode_scales)  # [3*n_valid]
        displacement = displacement_flat.reshape(n_valid, 3)  # [n_valid, 3]

        # Scale by pLDDT-based amplitude
        displacement *= amplitude_factors[:, None]

        # Apply displacement to all atoms in the affected residues
        new_coords = structure.coords.copy()
        for local_i, global_i in enumerate(valid_indices):
            new_coords[global_i] += displacement[local_i]

        conformations.append(new_coords.astype(np.float32))

    return conformations
