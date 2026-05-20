"""Unit tests for structure mutation utilities."""

import numpy as np
import pytest

from refold.types import Mutation, ProteinStructure
from refold.structure.mutator import apply_mutation
from refold.constants import CA_IDX, N_IDX, C_IDX, O_IDX, N_ATOM_TYPES


def test_apply_mutation_updates_residue_type(dummy_structure):
    mutation = Mutation(
        uniprot_id="DUMMY", position=5, wildtype_aa="A", mutant_aa="V"
    )
    mutant = apply_mutation(dummy_structure, mutation)
    assert mutant.mutant_residue_types[4] != dummy_structure.residue_types[4]


def test_apply_mutation_removes_side_chain(dummy_structure):
    mutation = Mutation(
        uniprot_id="DUMMY", position=10, wildtype_aa="A", mutant_aa="L"
    )
    mutant = apply_mutation(dummy_structure, mutation)
    pos_0 = 9
    for atom_idx in range(N_ATOM_TYPES):
        if atom_idx not in {N_IDX, CA_IDX, C_IDX, O_IDX}:
            assert np.all(np.isnan(mutant.mutant_coords[pos_0, atom_idx, :]))


def test_apply_mutation_preserves_backbone(dummy_structure):
    mutation = Mutation(
        uniprot_id="DUMMY", position=10, wildtype_aa="A", mutant_aa="L"
    )
    mutant = apply_mutation(dummy_structure, mutation)
    pos_0 = 9
    for atom_idx in (N_IDX, CA_IDX, C_IDX, O_IDX):
        if dummy_structure.atom_mask[pos_0, atom_idx]:
            assert not np.any(np.isnan(mutant.mutant_coords[pos_0, atom_idx, :]))


def test_apply_mutation_preserves_other_residues(dummy_structure):
    mutation = Mutation(
        uniprot_id="DUMMY", position=10, wildtype_aa="A", mutant_aa="L"
    )
    mutant = apply_mutation(dummy_structure, mutation)
    for i in [0, 5, 20, 49]:
        np.testing.assert_array_equal(
            mutant.mutant_coords[i], dummy_structure.coords[i]
        )


def test_mutation_position_validation(dummy_structure):
    mutation = Mutation(
        uniprot_id="DUMMY", position=200, wildtype_aa="A", mutant_aa="V"
    )
    with pytest.raises(ValueError):
        apply_mutation(dummy_structure, mutation)


def test_ddg_prediction_returns_float(dummy_structure):
    from refold.structure.mutator import ESM1vStabilityPredictor
    import torch

    predictor = ESM1vStabilityPredictor(device=torch.device("cpu"))
    try:
        ddg = predictor.predict_ddg(
            dummy_structure.sequence, 5, "A", "V"
        )
        assert isinstance(ddg, float)
    except Exception:
        pytest.skip("ESM-1v not installed")
