"""High-level predictor for single-mutation inference."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from refold.types import Mutation, REFOLDResult

logger = logging.getLogger(__name__)


class REFOLDPredictor:
    """
    Thin wrapper around REFOLDPipeline for programmatic use.

    Example::

        predictor = REFOLDPredictor.from_checkpoints("checkpoints/")
        result = predictor.predict("P06280", 152, "Y", "C")
    """

    def __init__(
        self,
        stage: int = 3,
        n_conformations: int = 50,
        n_molecules: int = 1000,
        device_str: Optional[str] = None,
    ):
        self.stage = stage
        self.n_conformations = n_conformations
        self.n_molecules = n_molecules
        self.device_str = device_str
        self._pipeline = None

    @classmethod
    def from_checkpoints(cls, checkpoint_dir: str | Path, **kwargs) -> "REFOLDPredictor":
        inst = cls(**kwargs)
        return inst

    def _get_pipeline(self):
        if self._pipeline is None:
            from refold.pipeline.single_mutation import REFOLDPipeline
            import torch

            device = None
            if self.device_str:
                device = torch.device(self.device_str)

            self._pipeline = REFOLDPipeline(
                device=device,
                n_conformations=self.n_conformations,
                n_molecules_per_pocket=self.n_molecules,
            )
        return self._pipeline

    def predict(
        self,
        uniprot_id: str,
        position: int,
        wildtype_aa: str,
        mutant_aa: str,
        gene_name: str = "",
        disease: str = "",
    ) -> REFOLDResult:
        mutation = Mutation(
            uniprot_id=uniprot_id,
            position=position,
            wildtype_aa=wildtype_aa,
            mutant_aa=mutant_aa,
            gene_name=gene_name,
            disease=disease,
        )
        return self._get_pipeline().run(mutation)

    def predict_batch(self, mutations: list[Mutation]) -> list[REFOLDResult]:
        pipeline = self._get_pipeline()
        results = []
        for mut in mutations:
            try:
                results.append(pipeline.run(mut))
            except Exception as e:
                logger.error(f"Prediction failed for {mut}: {e}")
        return results
