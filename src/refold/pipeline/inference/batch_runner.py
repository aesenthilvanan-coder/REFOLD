"""Batch runner for processing many mutations in parallel."""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import numpy as np

from refold.types import Mutation, REFOLDResult

logger = logging.getLogger(__name__)


class BatchRunner:
    """
    Run REFOLD pipeline on a list of mutations with optional parallelism.

    For GPU inference, set n_workers=1 to avoid GPU memory contention.
    For CPU/stage-1, n_workers=4 is reasonable.
    """

    def __init__(
        self,
        n_workers: int = 1,
        stage: int = 1,
        n_conformations: int = 10,
        n_molecules: int = 0,
        device_str: Optional[str] = None,
        output_dir: Optional[Path] = None,
    ):
        self.n_workers = n_workers
        self.stage = stage
        self.n_conformations = n_conformations
        self.n_molecules = n_molecules
        self.device_str = device_str
        self.output_dir = output_dir

    def run(
        self,
        mutations: list[Mutation],
        progress_callback=None,
    ) -> list[REFOLDResult]:
        from refold.pipeline.single_mutation import REFOLDPipeline
        try:
            import torch
            device = torch.device(self.device_str) if self.device_str else None
        except ImportError:
            device = None

        pipeline = REFOLDPipeline(
            device=device,
            n_conformations=self.n_conformations,
            n_molecules_per_pocket=self.n_molecules,
        )

        results = []
        t_start = time.time()

        if self.n_workers == 1:
            for i, mut in enumerate(mutations):
                try:
                    result = pipeline.run(mut)
                except Exception as e:
                    logger.error(f"Failed on {mut}: {e}")
                    result = _error_result(mut, str(e))
                results.append(result)
                if progress_callback:
                    progress_callback(i + 1, len(mutations), result)
                if self.output_dir:
                    _save_result(result, self.output_dir)
        else:
            # Multi-worker: each worker gets its own pipeline (no GPU sharing)
            def _run_one(mut: Mutation) -> REFOLDResult:
                from refold.pipeline.single_mutation import REFOLDPipeline
                p = REFOLDPipeline(
                    device=device,
                    n_conformations=self.n_conformations,
                    n_molecules_per_pocket=self.n_molecules,
                )
                try:
                    return p.run(mut)
                except Exception as e:
                    return _error_result(mut, str(e))

            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                futures = {executor.submit(_run_one, m): m for m in mutations}
                done_count = 0
                for future in as_completed(futures):
                    result = future.result()
                    results.append(result)
                    done_count += 1
                    if progress_callback:
                        progress_callback(done_count, len(mutations), result)
                    if self.output_dir:
                        _save_result(result, self.output_dir)

        elapsed = time.time() - t_start
        logger.info(f"Batch complete: {len(results)} mutations in {elapsed:.1f}s")
        return results


def _error_result(mutation: Mutation, error_msg: str) -> REFOLDResult:
    from refold.types import MutationClass, RescueAmenability
    return REFOLDResult(
        mutation=mutation,
        mutation_class=MutationClass.UNKNOWN,
        rescue_amenability=RescueAmenability.UNRESCUABLE,
        rescue_amenability_prob=0.0,
        rescue_probability=0.0,
        ddg_predicted=0.0,
        n_pockets_detected=0,
        n_molecules_generated=0,
        n_molecules_passing_filters=0,
        error_message=error_msg,
    )


def _save_result(result: REFOLDResult, output_dir: Path) -> None:
    import json
    from refold.utils.io import NumpyEncoder
    output_dir.mkdir(parents=True, exist_ok=True)
    key = result.mutation.to_key()
    path = output_dir / f"{key}.json"
    try:
        with open(path, "w") as f:
            json.dump(result.to_dict(), f, cls=NumpyEncoder, indent=2)
    except Exception as e:
        logger.warning(f"Failed to save result for {key}: {e}")
