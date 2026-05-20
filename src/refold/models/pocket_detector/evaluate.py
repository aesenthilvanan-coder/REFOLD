"""Evaluation metrics for the pocket detector."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def compute_pocket_detection_metrics(
    predicted_pockets: list,
    true_pocket_centers: list[np.ndarray],
    distance_threshold: float = 4.0,
) -> dict[str, float]:
    """
    Compute pocket detection metrics.

    A prediction is a true positive if any predicted pocket center is within
    *distance_threshold* Å of a true pocket center.

    Returns: precision, recall, f1
    """
    if not true_pocket_centers:
        return {"precision": 1.0, "recall": 1.0, "f1": 1.0, "n_pred": len(predicted_pockets)}

    if not predicted_pockets:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "n_pred": 0}

    pred_centers = np.array([p.center for p in predicted_pockets])

    tp = 0
    matched_preds = set()
    for true_center in true_pocket_centers:
        dists = np.sqrt(np.sum((pred_centers - true_center) ** 2, axis=1))
        best = int(np.argmin(dists))
        if dists[best] <= distance_threshold and best not in matched_preds:
            tp += 1
            matched_preds.add(best)

    precision = tp / len(predicted_pockets) if predicted_pockets else 0.0
    recall = tp / len(true_pocket_centers)
    f1 = 2 * precision * recall / (precision + recall + 1e-8)

    return {
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "n_pred": len(predicted_pockets),
        "n_true": len(true_pocket_centers),
        "tp": tp,
    }


def evaluate_druggability_ranking(
    pockets: list,
    true_druggable_ids: set[str],
) -> dict[str, float]:
    """Evaluate druggability ranking quality (AUROC proxy)."""
    if not pockets:
        return {"auroc": 0.5, "n_druggable": 0}

    scores = np.array([p.druggability_score for p in pockets])
    labels = np.array([1 if p.pocket_id in true_druggable_ids else 0 for p in pockets])

    if labels.sum() == 0 or labels.sum() == len(labels):
        return {"auroc": 0.5, "n_druggable": int(labels.sum())}

    # Compute AUROC via trapezoidal approximation
    order = np.argsort(-scores)
    labels_sorted = labels[order]
    n_pos = labels.sum()
    n_neg = len(labels) - n_pos

    tp, fp = 0, 0
    auc = 0.0
    prev_fp = 0

    for label in labels_sorted:
        if label == 1:
            tp += 1
        else:
            fp += 1
        if fp > prev_fp:
            auc += tp * (fp - prev_fp)
            prev_fp = fp

    auroc = auc / (n_pos * n_neg) if n_pos * n_neg > 0 else 0.5

    return {"auroc": float(np.clip(auroc, 0.0, 1.0)), "n_druggable": int(n_pos)}
