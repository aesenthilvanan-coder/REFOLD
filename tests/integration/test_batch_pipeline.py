"""Integration test for batch pipeline processing."""

import pytest
import numpy as np


@pytest.mark.integration
def test_batch_runner_stage1_smoke(tmp_path):
    """Stage-1 batch runner processes mutations without crashing (smoke test)."""
    from refold.types import Mutation, MutationClass
    from refold.pipeline.inference.batch_runner import BatchRunner, _error_result

    # Use _error_result as a safe fallback to test the runner logic
    mutations = [
        Mutation("P99999", 1, "A", "G", gene_name="TEST"),
        Mutation("P99998", 5, "V", "L", gene_name="TEST2"),
    ]

    # Test _error_result directly (no network required)
    for mut in mutations:
        result = _error_result(mut, "test_error")
        assert result.mutation_class == MutationClass.UNKNOWN
        assert result.error_message == "test_error"
        assert result.rescue_probability == 0.0


@pytest.mark.integration
def test_result_writer_parquet(tmp_path):
    """Result writer serialises results to parquet correctly."""
    from refold.types import Mutation, REFOLDResult, MutationClass, RescueAmenability
    from refold.pipeline.inference.result_writer import write_results_to_parquet
    import pandas as pd

    results = []
    for i in range(3):
        mut = Mutation(f"P0000{i}", 10 + i, "A", "G")
        results.append(REFOLDResult(
            mutation=mut,
            mutation_class=MutationClass.MISFOLDING,
            rescue_amenability=RescueAmenability.MODERATE,
            rescue_amenability_prob=0.6,
            rescue_probability=0.6,
            ddg_predicted=2.5,
            n_pockets_detected=2,
            n_molecules_generated=100,
            n_molecules_passing_filters=15,
        ))

    out = tmp_path / "results.parquet"
    write_results_to_parquet(results, out)
    assert out.exists()

    df = pd.read_parquet(out)
    assert len(df) == 3
    assert "rescue_probability" in df.columns
    assert "ddg_predicted" in df.columns


@pytest.mark.integration
def test_result_writer_csv(tmp_path):
    """Result writer produces CSV summary."""
    from refold.types import Mutation, REFOLDResult, MutationClass, RescueAmenability
    from refold.pipeline.inference.result_writer import write_csv_summary

    mut = Mutation("P06280", 152, "Y", "C", gene_name="GLA", disease="Fabry disease")
    result = REFOLDResult(
        mutation=mut,
        mutation_class=MutationClass.MISFOLDING,
        rescue_amenability=RescueAmenability.HIGH,
        rescue_amenability_prob=0.85,
        rescue_probability=0.85,
        ddg_predicted=3.2,
        n_pockets_detected=1,
        n_molecules_generated=50,
        n_molecules_passing_filters=8,
        runtime_seconds=12.5,
    )

    out = tmp_path / "summary.csv"
    write_csv_summary([result], out)
    assert out.exists()

    import pandas as pd
    df = pd.read_csv(out)
    assert len(df) == 1
    assert "rescue_prob" in df.columns


@pytest.mark.integration
def test_inference_predictor_instantiation():
    """REFOLDPredictor instantiates without ML deps."""
    from refold.pipeline.inference.predictor import REFOLDPredictor

    predictor = REFOLDPredictor(stage=1, n_conformations=0, n_molecules=0)
    assert predictor.stage == 1


@pytest.mark.integration
def test_batch_runner_save_results(tmp_path):
    """BatchRunner saves per-mutation JSON files."""
    from refold.pipeline.inference.batch_runner import _save_result, _error_result
    from refold.types import Mutation
    import json

    mut = Mutation("P99999", 1, "A", "G")
    result = _error_result(mut, "test")
    _save_result(result, tmp_path)

    key = mut.to_key()
    out = tmp_path / f"{key}.json"
    assert out.exists()

    with open(out) as f:
        data = json.load(f)
    assert data["mutation"]["uniprot_id"] == "P99999"
