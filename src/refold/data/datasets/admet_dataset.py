"""Dataset for ADMET predictor training."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.models.admet_predictor.evaluate import ADMET_TASK_NAMES

logger = logging.getLogger(__name__)


class ADMETDataset:
    """
    Dataset of (SMILES, labels) for multi-task ADMET training.

    Labels shape: [N, 8] float32, NaN for missing tasks.
    """

    def __init__(
        self,
        smiles_list: list[str],
        labels: np.ndarray,
        precompute_features: bool = False,
    ):
        self.smiles_list = smiles_list
        self.labels = labels
        self._features: Optional[list] = None

        assert len(smiles_list) == len(labels), "SMILES/labels length mismatch"
        assert labels.shape[1] == len(ADMET_TASK_NAMES), (
            f"Expected {len(ADMET_TASK_NAMES)} tasks, got {labels.shape[1]}"
        )

        if precompute_features:
            self._precompute()

    def _precompute(self) -> None:
        from refold.models.admet_predictor.features import build_admet_features
        logger.info(f"Precomputing ADMET features for {len(self.smiles_list)} molecules...")
        self._features = [build_admet_features(s) for s in self.smiles_list]

    def __len__(self) -> int:
        return len(self.smiles_list)

    def __getitem__(self, idx: int) -> dict:
        smiles = self.smiles_list[idx]
        label = self.labels[idx].copy()  # [8]

        if self._features is not None and self._features[idx] is not None:
            feat = self._features[idx]
        else:
            from refold.models.admet_predictor.features import build_admet_features
            feat = build_admet_features(smiles)

        if feat is None:
            from refold.models.admet_predictor.features import ADMET_FEAT_DIM
            feat = np.zeros(ADMET_FEAT_DIM, dtype=np.float32)

        return {
            "features": feat,
            "labels": label,
            "smiles": smiles,
        }

    def get_task_prevalences(self) -> dict[str, float]:
        """Return positive class fraction per task (ignoring NaN)."""
        result = {}
        for t, name in enumerate(ADMET_TASK_NAMES):
            col = self.labels[:, t]
            valid = col[~np.isnan(col)]
            if len(valid) > 0:
                result[name] = float(valid.mean())
            else:
                result[name] = 0.5
        return result
