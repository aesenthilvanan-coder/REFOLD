"""Shared pytest fixtures for REFOLD tests."""

import numpy as np
import pytest

from refold.types import Mutation, ProteinStructure, Pocket, PocketType


@pytest.fixture
def fab_mutation() -> Mutation:
    """GLA Y152C — Fabry disease chaperone benchmark variant."""
    return Mutation(
        uniprot_id="P06280",
        position=152,
        wildtype_aa="Y",
        mutant_aa="C",
        gene_name="GLA",
        disease="Fabry disease",
    )


@pytest.fixture
def tp53_mutation() -> Mutation:
    """TP53 R175H — Li-Fraumeni syndrome benchmark variant."""
    return Mutation(
        uniprot_id="P04637",
        position=175,
        wildtype_aa="R",
        mutant_aa="H",
        gene_name="TP53",
        disease="Li-Fraumeni syndrome",
    )


@pytest.fixture
def dummy_structure() -> ProteinStructure:
    """50-residue ProteinStructure with backbone atoms for testing."""
    n_res = 50
    n_atom_types = 37
    sequence = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM"[:n_res]

    coords = np.full((n_res, n_atom_types, 3), np.nan, dtype=np.float32)
    atom_mask = np.zeros((n_res, n_atom_types), dtype=bool)

    # Set backbone atoms: N=0, CA=1, C=2, O=3
    rng = np.random.default_rng(42)
    for i in range(n_res):
        base_pos = np.array([i * 3.8, 0.0, 0.0], dtype=np.float32)
        for atom_idx, offset in [(0, [0.0, 0.0, 0.0]), (1, [1.5, 0.0, 0.0]),
                                  (2, [2.5, 1.0, 0.0]), (3, [3.5, 0.0, 0.0])]:
            coords[i, atom_idx] = base_pos + np.array(offset, dtype=np.float32)
            atom_mask[i, atom_idx] = True

    residue_types = np.zeros(n_res, dtype=np.int64)
    residue_mask = atom_mask[:, 1]  # CA mask
    bfactors = rng.uniform(60.0, 95.0, n_res).astype(np.float32)

    return ProteinStructure(
        uniprot_id="DUMMY",
        sequence=sequence,
        coords=coords,
        residue_types=residue_types,
        residue_mask=residue_mask,
        atom_mask=atom_mask,
        bfactors=bfactors,
    )


@pytest.fixture
def dummy_pocket() -> Pocket:
    """TRANSIENT_MISFOLDING pocket for testing."""
    n_spheres = 30
    rng = np.random.default_rng(0)
    center = np.array([20.0, 5.0, 0.0], dtype=np.float32)
    alpha_spheres = (
        rng.normal(0, 3, (n_spheres, 3)).astype(np.float32) + center
    )

    return Pocket(
        pocket_id="pocket_001",
        pocket_type=PocketType.TRANSIENT_MISFOLDING,
        center=center,
        volume=650.0,
        druggability_score=0.75,
        hydrophobicity=0.65,
        residue_indices=list(range(5, 20)),
        alpha_sphere_coords=alpha_spheres,
        detection_frequency=0.65,
    )
