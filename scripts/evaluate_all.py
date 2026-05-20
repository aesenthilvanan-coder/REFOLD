#!/usr/bin/env python3
"""Evaluate REFOLD against benchmark targets."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.constants import BENCHMARK_MUTATIONS, RESULTS_DIR, PROCESSED_DIR

logger = logging.getLogger(__name__)

BENCHMARK_TARGETS = {
    "auroc": 0.85,
    "auprc": 0.70,
    "mcc": 0.55,
    "f1": 0.65,
    "ddg_pearson_r": 0.60,
}


def evaluate_on_benchmark(device_str: str | None = None) -> dict:
    """Run pipeline on all benchmark mutations and evaluate."""
    from refold.pipeline.single_mutation import REFOLDPipeline
    from refold.types import Mutation
    import torch

    device = torch.device(device_str) if device_str else None

    pipeline = REFOLDPipeline(
        device=device,
        n_conformations=10,
        n_molecules_per_pocket=0,
    )

    y_true, y_prob, ddg_pred, ddg_true = [], [], [], []

    for key in BENCHMARK_MUTATIONS:
        parts = key.split("_")
        if len(parts) != 4:
            continue
        uniprot_id, pos_str, wt, mt = parts
        try:
            mutation = Mutation(
                uniprot_id=uniprot_id,
                position=int(pos_str),
                wildtype_aa=wt,
                mutant_aa=mt,
            )
            result = pipeline.run(mutation)

            # Known chaperone-amenable mutations are positive
            y_true.append(1)
            y_prob.append(result.rescue_probability)
            ddg_pred.append(result.ddg_predicted)
        except Exception as e:
            logger.warning(f"Benchmark evaluation failed for {key}: {e}")

    if not y_true:
        logger.error("No benchmark results collected")
        return {}

    from refold.models.rescue_classifier.evaluate import compute_metrics
    metrics = compute_metrics(y_true, y_prob)

    if ddg_pred and ddg_true:
        import numpy as np
        from scipy.stats import pearsonr
        r, _ = pearsonr(ddg_pred, ddg_true)
        metrics["ddg_pearson_r"] = float(r)

    logger.info(f"Benchmark metrics: {metrics}")
    return metrics


def print_benchmark_report(metrics: dict) -> None:
    print("\n" + "=" * 60)
    print("REFOLD BENCHMARK EVALUATION REPORT")
    print("=" * 60)

    all_pass = True
    for metric, target in BENCHMARK_TARGETS.items():
        value = metrics.get(metric, 0.0)
        status = "PASS" if value >= target else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {metric:<20} {value:.3f}  (target: ≥{target:.2f})  [{status}]")

    print("=" * 60)
    print(f"Overall: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate REFOLD on benchmark")
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, default=RESULTS_DIR / "benchmark_report.json")
    args = parser.parse_args()

    setup_logging()
    metrics = evaluate_on_benchmark(args.device)

    if metrics:
        print_benchmark_report(metrics)
        import json
        from refold.utils.io import NumpyEncoder
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(metrics, f, indent=2, cls=NumpyEncoder)
        logger.info(f"Report saved to {args.output}")


if __name__ == "__main__":
    main()
