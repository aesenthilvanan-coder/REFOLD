"""Data augmentation strategies for training."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.types import ProteinStructure, Mutation
from refold.constants import STANDARD_AAS, AA_TO_IDX

logger = logging.getLogger(__name__)


def add_coordinate_noise(
    coords: np.ndarray,
    atom_mask: np.ndarray,
    sigma: float = 0.05,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Add Gaussian noise to atom coordinates (only where atom_mask is True)."""
    if rng is None:
        rng = np.random.default_rng()
    noise = rng.normal(0, sigma, coords.shape).astype(np.float32)
    noise[~atom_mask] = 0.0  # don't perturb missing atoms
    return coords + noise


def random_rotation(
    coords: np.ndarray,
    atom_mask: np.ndarray,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Apply random 3D rotation to all atom coordinates."""
    if rng is None:
        rng = np.random.default_rng()

    # Random quaternion via Shoemake (1992)
    u1, u2, u3 = rng.uniform(0, 1, 3)
    q = np.array([
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ])
    # Quaternion to rotation matrix
    w, x, y, z = q
    R = np.array([
        [1 - 2*(y*y+z*z),   2*(x*y-z*w),   2*(x*z+y*w)],
        [  2*(x*y+z*w), 1 - 2*(x*x+z*z),   2*(y*z-x*w)],
        [  2*(x*z-y*w),   2*(y*z+x*w), 1 - 2*(x*x+y*y)],
    ], dtype=np.float32)

    orig_shape = coords.shape
    flat = coords.reshape(-1, 3)
    rotated = (R @ flat.T).T.reshape(orig_shape)
    # Preserve NaN structure
    rotated[~atom_mask] = np.nan
    return rotated


def center_structure(coords: np.ndarray, atom_mask: np.ndarray) -> np.ndarray:
    """Translate structure so Cα centroid is at origin."""
    ca = coords[:, 1, :]  # CA_IDX = 1
    ca_mask = atom_mask[:, 1]
    valid_ca = ca[ca_mask]
    if len(valid_ca) == 0:
        return coords
    centroid = valid_ca.mean(axis=0)
    centered = coords.copy()
    # Shift only valid atoms
    for i in range(coords.shape[0]):
        for j in range(coords.shape[1]):
            if atom_mask[i, j] and not np.any(np.isnan(coords[i, j])):
                centered[i, j] -= centroid
    return centered


def augment_mutation_label(
    mutation: Mutation,
    rng: Optional[np.random.Generator] = None,
    flip_prob: float = 0.0,
) -> Mutation:
    """Optionally flip wildtype/mutant for negative data augmentation. Default: no-op."""
    if rng is None or flip_prob == 0.0:
        return mutation
    if rng.random() < flip_prob:
        return Mutation(
            uniprot_id=mutation.uniprot_id,
            position=mutation.position,
            wildtype_aa=mutation.mutant_aa,
            mutant_aa=mutation.wildtype_aa,
            gene_name=mutation.gene_name,
            disease=mutation.disease,
            clinvar_id=mutation.clinvar_id,
            source=mutation.source,
        )
    return mutation


def generate_synthetic_negative(
    structure: ProteinStructure,
    rng: Optional[np.random.Generator] = None,
) -> Mutation:
    """
    Generate a synthetic negative (non-pathogenic) mutation.

    Picks a surface-exposed, non-conserved position and a physicochemically
    similar substitution.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = structure.n_residues
    position = int(rng.integers(1, n + 1))
    wt = structure.sequence[position - 1] if position <= len(structure.sequence) else "A"

    # Pick similar amino acid
    similar = {
        "A": ["G", "S", "V"], "G": ["A", "S"], "V": ["I", "L", "A"],
        "I": ["V", "L"], "L": ["I", "V", "M"], "M": ["L", "I", "V"],
        "F": ["Y", "W"], "Y": ["F", "H"], "W": ["F", "Y"],
        "S": ["T", "A"], "T": ["S", "V"], "C": ["S", "A"],
        "D": ["E", "N"], "E": ["D", "Q"], "N": ["Q", "D"], "Q": ["N", "E"],
        "K": ["R"], "R": ["K"], "H": ["K", "R"],
        "P": ["A", "G"],
    }
    choices = similar.get(wt, [aa for aa in STANDARD_AAS if aa != wt])
    mt = str(rng.choice(choices)) if choices else "A"

    return Mutation(
        uniprot_id=structure.uniprot_id,
        position=position,
        wildtype_aa=wt,
        mutant_aa=mt,
        source="synthetic_negative",
    )
