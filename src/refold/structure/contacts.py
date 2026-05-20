"""Residue contact detection from Cα distances."""

from __future__ import annotations

import numpy as np

from refold.constants import CONTACT_CA_CUTOFF, INTERFACE_HEAVY_ATOM_CUTOFF
from refold.types import ProteinStructure


def compute_ca_distance_matrix(structure: ProteinStructure) -> np.ndarray:
    """Return [N_res, N_res] float32 pairwise Cα distances; NaN for missing residues."""
    ca = structure.ca_coords  # [N, 3]
    diff = ca[:, None, :] - ca[None, :, :]  # [N, N, 3]
    dist = np.sqrt(np.sum(diff ** 2, axis=-1)).astype(np.float32)  # [N, N]
    # Mask out unresolved residues
    missing = np.any(np.isnan(ca), axis=-1)  # [N]
    dist[missing, :] = np.nan
    dist[:, missing] = np.nan
    return dist


def compute_contact_map(
    structure: ProteinStructure,
    cutoff: float = CONTACT_CA_CUTOFF,
    min_seq_sep: int = 6,
) -> np.ndarray:
    """Return [N_res, N_res] bool contact map."""
    dist = compute_ca_distance_matrix(structure)
    contacts = (dist <= cutoff) & ~np.isnan(dist)
    # Exclude contacts between sequence-adjacent residues
    n = len(structure.sequence)
    for k in range(-min_seq_sep + 1, min_seq_sep):
        idx = np.arange(max(0, k), min(n, n + k))
        jdx = idx - k
        valid = (idx >= 0) & (idx < n) & (jdx >= 0) & (jdx < n)
        contacts[idx[valid], jdx[valid]] = False
    return contacts


def get_residues_within_radius(
    structure: ProteinStructure,
    center_idx: int,
    radius: float = CONTACT_CA_CUTOFF,
) -> list[int]:
    """Return 0-based indices of residues with Cα within *radius* Å of residue *center_idx*."""
    dist = compute_ca_distance_matrix(structure)
    row = dist[center_idx]
    return [int(i) for i in np.where((row <= radius) & ~np.isnan(row))[0] if i != center_idx]


def compute_interface_contacts(
    structure_a: ProteinStructure,
    structure_b: ProteinStructure,
    cutoff: float = INTERFACE_HEAVY_ATOM_CUTOFF,
) -> tuple[list[int], list[int]]:
    """Return (residue_indices_a, residue_indices_b) at the interface of two structures."""
    ca_a = structure_a.ca_coords  # [Na, 3]
    ca_b = structure_b.ca_coords  # [Nb, 3]

    diff = ca_a[:, None, :] - ca_b[None, :, :]  # [Na, Nb, 3]
    dist = np.sqrt(np.sum(diff ** 2, axis=-1))   # [Na, Nb]

    # Mask missing
    missing_a = np.any(np.isnan(ca_a), axis=-1)
    missing_b = np.any(np.isnan(ca_b), axis=-1)
    dist[missing_a, :] = np.inf
    dist[:, missing_b] = np.inf

    contact_mask = dist <= cutoff
    res_a = sorted(set(int(i) for i in np.where(contact_mask.any(axis=1))[0]))
    res_b = sorted(set(int(j) for j in np.where(contact_mask.any(axis=0))[0]))
    return res_a, res_b
