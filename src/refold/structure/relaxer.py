"""Lightweight structure relaxation via gradient descent on clashes."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.constants import CA_IDX, N_IDX, C_IDX, O_IDX

logger = logging.getLogger(__name__)

# Ideal bond lengths (Å)
_BOND_LENGTHS = {
    ("N", "CA"): 1.46,
    ("CA", "C"): 1.52,
    ("C", "N"): 1.33,   # peptide bond
    ("C", "O"): 1.23,
}

# Minimum allowed Cα–Cα distance for non-adjacent residues (Å)
_MIN_CA_CA_DIST = 3.8


def relax_backbone(
    coords: np.ndarray,
    atom_mask: np.ndarray,
    n_steps: int = 50,
    step_size: float = 0.01,
) -> np.ndarray:
    """
    Minimal gradient-descent relaxation to resolve steric clashes in backbone.

    Operates only on N, CA, C, O atoms (indices 0-4).

    Args:
        coords: [N_res, 37, 3] float32
        atom_mask: [N_res, 37] bool
        n_steps: number of gradient steps
        step_size: learning rate

    Returns:
        [N_res, 37, 3] relaxed coords (same dtype)
    """
    coords = coords.copy().astype(np.float64)
    n_res = coords.shape[0]

    for _ in range(n_steps):
        grad = np.zeros_like(coords)

        # Clash repulsion between Cα of non-adjacent residues
        ca = coords[:, CA_IDX, :]  # [N, 3]
        for i in range(n_res):
            if not atom_mask[i, CA_IDX]:
                continue
            for j in range(i + 2, min(i + 20, n_res)):
                if not atom_mask[j, CA_IDX]:
                    continue
                diff = ca[i] - ca[j]
                dist = float(np.linalg.norm(diff))
                if dist < _MIN_CA_CA_DIST and dist > 1e-6:
                    force = ((_MIN_CA_CA_DIST - dist) / dist) * diff
                    grad[i, CA_IDX] += force
                    grad[j, CA_IDX] -= force

        coords -= step_size * grad

    return coords.astype(np.float32)


def check_clashes(coords: np.ndarray, atom_mask: np.ndarray) -> int:
    """Return number of severe Cα–Cα clashes (< 3.0 Å) between non-adjacent residues."""
    ca = coords[:, CA_IDX, :]
    n = len(ca)
    n_clashes = 0
    for i in range(n):
        if not atom_mask[i, CA_IDX]:
            continue
        for j in range(i + 2, n):
            if not atom_mask[j, CA_IDX]:
                continue
            dist = float(np.linalg.norm(ca[i] - ca[j]))
            if dist < 3.0:
                n_clashes += 1
    return n_clashes
