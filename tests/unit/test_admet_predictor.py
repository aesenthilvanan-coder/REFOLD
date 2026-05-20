"""Unit tests for ADMET predictor."""

import numpy as np
import pytest


def test_build_admet_features_valid_smiles():
    pytest.importorskip("rdkit")
    from refold.models.admet_predictor.features import build_admet_features, ADMET_FEAT_DIM

    feat = build_admet_features("CCO")
    assert feat is not None
    assert feat.shape == (ADMET_FEAT_DIM,)
    assert feat.dtype == np.float32
    assert np.all(feat >= 0)


def test_build_admet_features_invalid_smiles():
    pytest.importorskip("rdkit")
    from refold.models.admet_predictor.features import build_admet_features

    feat = build_admet_features("not_a_smiles_XXXXXXX")
    assert feat is None


def test_morgan_fp_dimensions():
    pytest.importorskip("rdkit")
    from refold.models.admet_predictor.features import smiles_to_morgan_fp

    fp = smiles_to_morgan_fp("c1ccccc1")
    assert fp is not None
    assert fp.shape == (2048,)
    assert set(fp).issubset({0.0, 1.0})


def test_maccs_keys_dimensions():
    pytest.importorskip("rdkit")
    from refold.models.admet_predictor.features import smiles_to_maccs_keys

    fp = smiles_to_maccs_keys("c1ccccc1")
    assert fp is not None
    assert fp.shape == (167,)


def test_admet_evaluate_auroc():
    from refold.models.admet_predictor.evaluate import _compute_auroc

    y_true = np.array([0, 0, 1, 1], dtype=float)
    y_score = np.array([0.1, 0.2, 0.8, 0.9])
    auc = _compute_auroc(y_true, y_score)
    assert abs(auc - 1.0) < 0.01


def test_admet_evaluate_random_auroc():
    from refold.models.admet_predictor.evaluate import _compute_auroc

    rng = np.random.default_rng(42)
    y_true = rng.integers(0, 2, 100).astype(float)
    y_score = rng.uniform(0, 1, 100)
    auc = _compute_auroc(y_true, y_score)
    assert 0.0 <= auc <= 1.0


def test_admet_compute_metrics_shape():
    from refold.models.admet_predictor.evaluate import compute_admet_metrics, ADMET_TASK_NAMES

    n = 50
    n_tasks = len(ADMET_TASK_NAMES)
    rng = np.random.default_rng(0)
    y_true = rng.integers(0, 2, (n, n_tasks)).astype(float)
    y_pred = rng.uniform(0, 1, (n, n_tasks)).astype(float)

    results = compute_admet_metrics(y_true, y_pred)
    assert len(results) == n_tasks
    for task, metrics in results.items():
        assert "auroc" in metrics
        assert 0.0 <= metrics["auroc"] <= 1.0


def test_synthesizability_score_range():
    pytest.importorskip("rdkit")
    from refold.scoring.synthesizability import compute_synthetic_accessibility

    smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    for smi in smiles_list:
        score = compute_synthetic_accessibility(smi)
        assert 0.0 <= score <= 1.0


def test_druggability_score_range():
    from refold.scoring.druggability import compute_druggability_score
    from refold.types import PocketType

    score = compute_druggability_score(
        volume=600.0,
        fpocket_score=0.7,
        hydrophobicity=0.6,
        detection_freq=0.8,
        pocket_type=PocketType.TRANSIENT_MISFOLDING,
        n_residues=15,
    )
    assert 0.0 <= score <= 1.0


def test_binding_affinity_is_negative():
    pytest.importorskip("rdkit")
    from refold.scoring.binding_affinity import score_molecule_pocket_affinity
    from refold.types import Pocket, PocketType
    import numpy as np

    pocket = Pocket(
        pocket_id="P1",
        pocket_type=PocketType.TRANSIENT_MISFOLDING,
        center=np.zeros(3),
        volume=600.0,
        druggability_score=0.7,
        hydrophobicity=0.6,
        residue_indices=[0, 1, 2],
        alpha_sphere_coords=np.zeros((3, 3)),
        detection_frequency=0.8,
    )
    affinity = score_molecule_pocket_affinity("CCO", pocket)
    assert affinity <= 0.0
    assert -15.0 <= affinity <= 0.0
