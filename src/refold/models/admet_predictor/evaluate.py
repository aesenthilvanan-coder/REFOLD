"""Evaluation metrics for the ADMET predictor."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

ADMET_TASK_NAMES = [
    "herg_inhibition",
    "bbb_permeability",
    "cyp3a4_inhibition",
    "oral_bioavailability",
    "solubility",
    "hia_absorption",
    "pgp_substrate",
    "ames_mutagenicity",
]


def compute_admet_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    task_names: Optional[list[str]] = None,
) -> dict[str, dict[str, float]]:
    """
    Compute per-task AUROC and AUPRC.

    Args:
        y_true: [N, n_tasks] float32 (0/1, NaN for missing)
        y_pred: [N, n_tasks] float32 (probabilities)

    Returns:
        dict task_name → {"auroc": float, "auprc": float, "n_samples": int}
    """
    if task_names is None:
        task_names = ADMET_TASK_NAMES

    n_tasks = y_true.shape[1] if y_true.ndim > 1 else 1
    results = {}

    for t in range(min(n_tasks, len(task_names))):
        task = task_names[t]
        yt = y_true[:, t] if y_true.ndim > 1 else y_true
        yp = y_pred[:, t] if y_pred.ndim > 1 else y_pred

        # Remove NaN
        mask = ~np.isnan(yt)
        yt_clean = yt[mask]
        yp_clean = yp[mask]

        if len(yt_clean) < 2 or yt_clean.sum() == 0:
            results[task] = {"auroc": 0.5, "auprc": 0.0, "n_samples": int(mask.sum())}
            continue

        auroc = _compute_auroc(yt_clean, yp_clean)
        auprc = _compute_auprc(yt_clean, yp_clean)

        results[task] = {
            "auroc": auroc,
            "auprc": auprc,
            "n_samples": int(mask.sum()),
        }

    return results


def _compute_auroc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUROC via sorted pairs."""
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    n_pos = y_sorted.sum()
    n_neg = len(y_sorted) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    tp_cum = np.cumsum(y_sorted)
    fp_cum = np.cumsum(1 - y_sorted)
    tpr = tp_cum / n_pos
    fpr = fp_cum / n_neg
    auc = float(np.trapezoid(tpr, fpr))
    return float(np.clip(auc, 0.0, 1.0))


def _compute_auprc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute area under the precision-recall curve."""
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    n_pos = y_sorted.sum()
    if n_pos == 0:
        return 0.0
    tp_cum = np.cumsum(y_sorted)
    precision = tp_cum / (np.arange(len(y_sorted)) + 1)
    recall = tp_cum / n_pos
    auc = float(np.trapezoid(precision, recall))
    return float(np.clip(auc, 0.0, 1.0))
