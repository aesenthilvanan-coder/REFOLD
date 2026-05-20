"""Training loop for the multi-task ADMET predictor."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from refold.constants import CHECKPOINT_DIR, ADMET_N_TASKS

logger = logging.getLogger(__name__)


def train_admet_predictor(
    train_smiles: list[str],
    train_labels: np.ndarray,
    val_smiles: list[str],
    val_labels: np.ndarray,
    n_epochs: int = 100,
    lr: float = 1e-4,
    batch_size: int = 256,
    checkpoint_dir: Path = CHECKPOINT_DIR / "admet_predictor",
    device_str: Optional[str] = None,
) -> dict:
    """
    Train the multi-task ADMET predictor.

    Args:
        train_labels: [N_train, 8] float32 — one column per ADMET task
        val_labels: [N_val, 8] float32
    """
    try:
        import torch
        from torch.optim import AdamW
        from torch.optim.lr_scheduler import CosineAnnealingLR
        from refold.models.admet_predictor.model import ADMETPredictor
        from refold.models.admet_predictor.features import build_admet_features
        from refold.utils.device import get_device
    except ImportError as e:
        raise ImportError(f"ML dependencies required: {e}")

    device = get_device(device_str)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    # Pre-compute features
    logger.info("Computing molecular fingerprints for training set...")
    train_feats = [build_admet_features(s) for s in train_smiles]
    val_feats = [build_admet_features(s) for s in val_smiles]

    # Filter out failed
    train_valid = [(f, l) for f, l in zip(train_feats, train_labels) if f is not None]
    val_valid = [(f, l) for f, l in zip(val_feats, val_labels) if f is not None]

    if not train_valid:
        logger.error("No valid training molecules")
        return {}

    train_x = torch.tensor(np.array([f for f, _ in train_valid]), dtype=torch.float32)
    train_y = torch.tensor(np.array([l for _, l in train_valid]), dtype=torch.float32)
    val_x = torch.tensor(np.array([f for f, _ in val_valid]), dtype=torch.float32)
    val_y = torch.tensor(np.array([l for _, l in val_valid]), dtype=torch.float32)

    from refold.models.admet_predictor.features import ADMET_FEAT_DIM
    model = ADMETPredictor(input_dim=ADMET_FEAT_DIM, n_tasks=ADMET_N_TASKS).to(device)
    optimizer = AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    n_steps = n_epochs * max(1, len(train_valid) // batch_size)
    scheduler = CosineAnnealingLR(optimizer, T_max=n_steps)

    criterion = torch.nn.BCEWithLogitsLoss()
    best_val_loss = float("inf")
    metrics = {"train_losses": [], "val_losses": [], "best_epoch": 0}

    for epoch in range(n_epochs):
        model.train()
        perm = torch.randperm(len(train_x))
        epoch_losses = []

        for start in range(0, len(train_x), batch_size):
            idx = perm[start:start + batch_size]
            xb = train_x[idx].to(device)
            yb = train_y[idx].to(device)

            optimizer.zero_grad()
            logits = model(xb)
            # Mask NaN labels
            mask = ~torch.isnan(yb)
            loss = criterion(logits[mask], yb[mask])
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            scheduler.step()
            epoch_losses.append(float(loss.item()))

        # Validation
        model.eval()
        with torch.no_grad():
            val_logits = model(val_x.to(device))
            mask = ~torch.isnan(val_y.to(device))
            val_loss = float(criterion(val_logits[mask], val_y.to(device)[mask]).item())

        train_loss = float(np.mean(epoch_losses)) if epoch_losses else float("inf")
        metrics["train_losses"].append(train_loss)
        metrics["val_losses"].append(val_loss)

        logger.info(f"Epoch {epoch + 1}/{n_epochs}  train={train_loss:.4f}  val={val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            metrics["best_epoch"] = epoch + 1
            torch.save(
                {"epoch": epoch, "model_state": model.state_dict(), "val_loss": val_loss},
                checkpoint_dir / "best.pt",
            )

        torch.save(
            {"epoch": epoch, "model_state": model.state_dict()},
            checkpoint_dir / "latest.pt",
        )

    return metrics
