"""ANM/ENM conformational ensemble generation.

Generates 50 partially unfolded states per mutation using
Anisotropic Network Model normal mode analysis.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from refold.constants import (
    CA_IDX, ENM_CUTOFF, N_ENM_MODES, N_CONFORMATIONS,
    ENM_AMPLITUDE_SCALE, ENM_N_SKIP_TRIVIAL,
)

if TYPE_CHECKING:
    from refold.types import ProteinStructure, MutantStructure

logger = logging.getLogger(__name__)


def compute_anm_hessian(
    ca_coords: np.ndarray,
    cutoff: float = ENM_CUTOFF,
    spring_constant: float = 1.0,
) -> np.ndarray:
    """Compute the 3N×3N ANM Hessian matrix.

    Args:
        ca_coords: [N_res, 3] Cα coordinates.
        cutoff: Distance cutoff for ANM contacts (Å).
        spring_constant: Spring constant γ.

    Returns:
        Hessian matrix [3N, 3N].
    """
    n_res = len(ca_coords)
    H = np.zeros((3 * n_res, 3 * n_res), dtype=np.float64)

    for i in range(n_res):
        for j in range(i + 1, n_res):
            r_ij = ca_coords[j] - ca_coords[i]
            dist = np.linalg.norm(r_ij)
            if dist > cutoff or dist < 1e-6:
                continue

            r_unit = r_ij / dist
            # Outer product of unit vector
            H_ij = -spring_constant * np.outer(r_unit, r_unit)

            # Off-diagonal blocks
            H[3*i:3*i+3, 3*j:3*j+3] = H_ij
            H[3*j:3*j+3, 3*i:3*i+3] = H_ij

            # Diagonal blocks
            H[3*i:3*i+3, 3*i:3*i+3] -= H_ij
            H[3*j:3*j+3, 3*j:3*j+3] -= H_ij

    return H


def generate_conformational_ensemble(
    structure: "ProteinStructure",
    n_conformations: int = N_CONFORMATIONS,
    n_modes: int = N_ENM_MODES,
    amplitude_scale: float = ENM_AMPLITUDE_SCALE,
    random_seed: int = 42,
) -> list[np.ndarray]:
    """Generate conformational ensemble via ANM normal modes.

    Eigendecomposes Hessian, skips 6 trivial modes (3 translation + 3 rotation),
    amplitude ∝ 1/√eigenvalue, pLDDT-based scaling.
    Conformation[0] is the original structure.

    Returns list of n_conformations Cα coordinate arrays [N_res, 3].
    """
    ca_coords = structure.ca_coords.copy()
    valid_mask = structure.residue_mask & ~np.any(np.isnan(ca_coords), axis=-1)
    ca_valid = ca_coords[valid_mask]
    n_res = len(ca_valid)

    if n_res < 5:
        logger.warning(f"Too few residues ({n_res}) for ENM — returning original only")
        return [ca_coords] * n_conformations

    try:
        H = compute_anm_hessian(ca_valid, cutoff=ENM_CUTOFF)
        eigenvalues, eigenvectors = np.linalg.eigh(H)
    except np.linalg.LinAlgError as e:
        logger.warning(f"ANM diagonalization failed: {e}")
        return [ca_coords] * n_conformations

    # Skip 6 trivial modes (near-zero eigenvalues: translation + rotation)
    skip = ENM_N_SKIP_TRIVIAL
    n_available = len(eigenvalues) - skip
    n_modes_use = min(n_modes, n_available)

    if n_modes_use <= 0:
        return [ca_coords] * n_conformations

    non_trivial_evals = eigenvalues[skip:skip + n_modes_use]
    non_trivial_evecs = eigenvectors[:, skip:skip + n_modes_use]

    # Amplitude ∝ 1/√eigenvalue
    amplitudes = amplitude_scale / np.sqrt(np.maximum(non_trivial_evals, 1e-10))

    # pLDDT-based scaling: low pLDDT residues get larger perturbations
    plddt = structure.bfactors[valid_mask]
    plddt_scale = np.clip(1.0 - plddt / 100.0, 0.1, 1.0)
    plddt_scale_3d = np.repeat(plddt_scale, 3)

    rng = np.random.default_rng(random_seed)
    conformations = [ca_coords.copy()]

    for _ in range(n_conformations - 1):
        coeffs = rng.normal(0, 1, n_modes_use) * amplitudes
        displacement_3d = (non_trivial_evecs * coeffs[None, :]).sum(axis=1)
        displacement_3d = displacement_3d * plddt_scale_3d

        new_ca = ca_valid.copy()
        new_ca += displacement_3d.reshape(n_res, 3)

        # Write back into full-size array
        new_ca_full = ca_coords.copy()
        new_ca_full[valid_mask] = new_ca
        conformations.append(new_ca_full)

    return conformations


def generate_mutant_ensemble(
    mutant_structure: "MutantStructure",
    n_conformations: int = N_CONFORMATIONS,
    n_modes: int = N_ENM_MODES,
    amplitude_scale: float = ENM_AMPLITUDE_SCALE,
    random_seed: int = 42,
) -> list[np.ndarray]:
    """Generate conformational ensemble for mutant — extra perturbation ∝ ΔΔG."""
    wt = mutant_structure.wildtype
    ddg = mutant_structure.ddg_predicted or 0.0

    # Extra perturbation proportional to destabilization
    ddg_scale = 1.0 + min(abs(ddg) / 5.0, 1.0)

    # Use mutant coordinates (side chains erased) as reference
    from refold.types import ProteinStructure
    mutant_as_structure = ProteinStructure(
        uniprot_id=wt.uniprot_id + "_mut",
        sequence=wt.sequence,
        coords=mutant_structure.mutant_coords,
        residue_types=mutant_structure.mutant_residue_types,
        residue_mask=wt.residue_mask,
        atom_mask=wt.atom_mask,
        bfactors=wt.bfactors,
        sse_ids=wt.sse_ids,
        phi_psi=wt.phi_psi,
        rel_asa=wt.rel_asa,
    )

    return generate_conformational_ensemble(
        mutant_as_structure,
        n_conformations=n_conformations,
        n_modes=n_modes,
        amplitude_scale=amplitude_scale * ddg_scale,
        random_seed=random_seed,
    )
