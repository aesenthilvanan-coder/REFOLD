"""Integration test for the single mutation pipeline."""

import pytest
from pathlib import Path
from refold.types import Mutation


@pytest.mark.integration
def test_fabry_pipeline_smoke(tmp_path):
    """Smoke test: run the GLA Y152C Fabry mutation through stage 1."""
    from refold.pipeline.single_mutation import REFOLDPipeline
    import torch

    mutation = Mutation(
        uniprot_id="P06280",
        position=152,
        wildtype_aa="Y",
        mutant_aa="C",
        gene_name="GLA",
        disease="Fabry disease",
    )

    pipeline = REFOLDPipeline(
        device=torch.device("cpu"),
        n_conformations=3,
        n_molecules_per_pocket=0,
    )

    result = pipeline.run(mutation, output_dir=tmp_path)
    assert result.mutation.uniprot_id == "P06280"
    assert result.ddg_predicted is not None
