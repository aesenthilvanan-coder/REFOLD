"""Evaluation metrics for the rescue classifier."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


def compute_metrics(
    y_true: list[int],
    y_prob: list[float],
    threshold: float = 0.5,
) -> dict[str, float]:
    """Compute AUROC, AUPRC, MCC, F1 for binary classification."""
    try:
        from sklearn.metrics import (
            roc_auc_score,
            average_precision_score,
            matthews_corrcoef,
            f1_score,
        )

        y_pred = [int(p >= threshold) for p in y_prob]
        return {
            "auroc": roc_auc_score(y_true, y_prob),
            "auprc": average_precision_score(y_true, y_prob),
            "mcc": matthews_corrcoef(y_true, y_pred),
            "f1": f1_score(y_true, y_pred, zero_division=0),
        }
    except Exception as e:
        logger.warning(f"Metrics computation failed: {e}")
        return {"auroc": 0.0, "auprc": 0.0, "mcc": 0.0, "f1": 0.0}


BENCHMARK_TARGETS = {
    "auroc": 0.85,
    "auprc": 0.70,
    "mcc": 0.55,
    "f1": 0.65,
}


def check_benchmark_targets(metrics: dict[str, float]) -> dict[str, bool]:
    return {
        metric: metrics.get(metric, 0.0) >= target
        for metric, target in BENCHMARK_TARGETS.items()
    }
