"""Unit tests for scoring and filtering utilities."""

import pytest

from refold.scoring.filters import (
    compute_lipinski_properties,
    compute_sa_score,
    compute_qed,
    check_pains,
)
from refold.scoring.rescue_probability import compute_final_rescue_probability
from refold.types import GeneratedMolecule


ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
INVALID_SMILES = "not_a_smiles"


def test_lipinski_known_drug_passes():
    props = compute_lipinski_properties(ASPIRIN_SMILES)
    assert props["passes_lipinski"] is True
    assert props["mw"] < 200.0
    assert props["hbd"] <= 5


def test_lipinski_invalid_smiles():
    props = compute_lipinski_properties(INVALID_SMILES)
    assert props["passes_lipinski"] is False
    assert props["mw"] == 0.0


def test_sa_score_aspirin():
    score = compute_sa_score(ASPIRIN_SMILES)
    assert 1.0 <= score <= 4.0


def test_qed_aspirin():
    qed = compute_qed(ASPIRIN_SMILES)
    assert 0.3 <= qed <= 1.0


def test_pains_aspirin_is_not_pains():
    is_pains = check_pains(ASPIRIN_SMILES)
    assert is_pains is False


def test_rescue_probability_range(dummy_pocket):
    mol = GeneratedMolecule(
        smiles=ASPIRIN_SMILES,
        pocket_id="pocket_001",
        predicted_affinity_kcal=-6.0,
        qed_score=0.7,
    )
    prob = compute_final_rescue_probability(
        mol, dummy_pocket, rescue_amenability_prob=0.8, ddg_protein=3.0
    )
    assert 0.0 <= prob <= 1.0


def test_high_affinity_increases_rescue_probability(dummy_pocket):
    mol_low = GeneratedMolecule(
        smiles=ASPIRIN_SMILES,
        pocket_id="pocket_001",
        predicted_affinity_kcal=-1.0,
        qed_score=0.5,
    )
    mol_high = GeneratedMolecule(
        smiles=ASPIRIN_SMILES,
        pocket_id="pocket_001",
        predicted_affinity_kcal=-10.0,
        qed_score=0.5,
    )
    prob_low = compute_final_rescue_probability(
        mol_low, dummy_pocket, rescue_amenability_prob=0.7, ddg_protein=3.0
    )
    prob_high = compute_final_rescue_probability(
        mol_high, dummy_pocket, rescue_amenability_prob=0.7, ddg_protein=3.0
    )
    assert prob_high > prob_low
