"""Training loop for the pocket detector GNN."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from refold.constants import CHECKPOINT_DIR

logger = logging.getLogger(__name__)


def train_pocket_detector(
    train_data: list,
    val_data: list,
    n_epochs: int = 50,
    lr: float = 1e-4,
    batch_size: int = 16,
    checkpoint_dir: Path = CHECKPOINT_DIR / "pocket_detector",
    device_str: Optional[str] = None,
) -> dict:
    """Train the pocket detector GNN model."""
    try:
        import torch
        from torch.optim import AdamW
        from torch.optim.lr_scheduler import CosineAnnealingLR
        from refold.models.pocket_detector.model import PocketDetectorGNN
        from refold.utils.device import get_device
    except ImportError as e:
        raise ImportError(f"ML dependencies required for training: {e}")

    device = get_device(device_str)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    model = PocketDetectorGNN().to(device)
    optimizer = AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = CosineAnnealingLR(optimizer, T_max=n_epochs * max(1, len(train_data) // batch_size))

    best_val_loss = float("inf")
    metrics = {"train_losses": [], "val_losses": [], "best_epoch": 0}

    for epoch in range(n_epochs):
        model.train()
        train_losses = []

        # Mini-batch training
        rng = np.random.default_rng(epoch)
        indices = rng.permutation(len(train_data))

        for start in range(0, len(indices), batch_size):
            batch_idx = indices[start:start + batch_size]
            batch = [train_data[i] for i in batch_idx]

            optimizer.zero_grad()
            loss = model.compute_loss(batch, device)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            scheduler.step()
            train_losses.append(float(loss.item()))

        # Validation
        model.eval()
        val_losses = []
        with torch.no_grad():
            for start in range(0, len(val_data), batch_size):
                batch = val_data[start:start + batch_size]
                loss = model.compute_loss(batch, device)
                val_losses.append(float(loss.item()))

        train_loss = np.mean(train_losses) if train_losses else float("inf")
        val_loss = np.mean(val_losses) if val_losses else float("inf")

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
            {"epoch": epoch, "model_state": model.state_dict(), "val_loss": val_loss},
            checkpoint_dir / "latest.pt",
        )

    return metrics
