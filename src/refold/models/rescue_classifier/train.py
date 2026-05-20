"""Rescue classifier training loop."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn

from refold.constants import (
    CHECKPOINT_DIR, EMA_DECAY, GRAD_CLIP_NORM,
)
from refold.models.rescue_classifier.model import RescueClassifier, FocalLoss

logger = logging.getLogger(__name__)


class EMA:
    """Exponential Moving Average for model parameters."""

    def __init__(self, model: nn.Module, decay: float = EMA_DECAY):
        self.decay = decay
        self.shadow: dict[str, torch.Tensor] = {}
        for name, param in model.named_parameters():
            if param.requires_grad:
                self.shadow[name] = param.data.clone()

    def update(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if param.requires_grad and name in self.shadow:
                self.shadow[name] = (
                    self.decay * self.shadow[name]
                    + (1.0 - self.decay) * param.data
                )

    def apply_shadow(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if name in self.shadow:
                param.data.copy_(self.shadow[name])

    def restore(self, model: nn.Module, original: dict[str, torch.Tensor]) -> None:
        for name, param in model.named_parameters():
            if name in original:
                param.data.copy_(original[name])


class RescueClassifierTrainer:
    """Training loop for the rescue classifier."""

    def __init__(
        self,
        model: RescueClassifier,
        device: torch.device,
        checkpoint_dir: Path | None = None,
        learning_rate: float = 3e-4,
        weight_decay: float = 0.01,
        gradient_accumulation: int = 4,
    ):
        self.model = model.to(device)
        self.device = device
        self.checkpoint_dir = checkpoint_dir or CHECKPOINT_DIR / "rescue_classifier"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=learning_rate, weight_decay=weight_decay
        )
        self.criterion = FocalLoss()
        self.ema = EMA(model)
        self.gradient_accumulation = gradient_accumulation
        self.best_val_loss = float("inf")
        self.global_step = 0

    def train_epoch(
        self,
        dataloader: Any,
        epoch: int,
    ) -> float:
        self.model.train()
        total_loss = 0.0
        n_batches = 0
        self.optimizer.zero_grad()

        for step, batch in enumerate(dataloader):
            batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                     for k, v in batch.items()}

            logits = self.model(
                batch["node_feats"],
                batch["edge_index"],
                batch["edge_feats"],
                batch["esm2_embedding"],
                batch["thermo_features"],
                batch["evo_features"],
                batch.get("node_mask"),
            )

            labels = batch["labels"].unsqueeze(-1)
            loss = self.criterion(logits, labels) / self.gradient_accumulation
            loss.backward()
            total_loss += loss.item() * self.gradient_accumulation

            if (step + 1) % self.gradient_accumulation == 0:
                nn.utils.clip_grad_norm_(self.model.parameters(), GRAD_CLIP_NORM)
                self.optimizer.step()
                self.optimizer.zero_grad()
                self.ema.update(self.model)
                self.global_step += 1

            n_batches += 1

        return total_loss / max(n_batches, 1)

    @torch.no_grad()
    def evaluate(self, dataloader: Any) -> dict[str, float]:
        self.model.eval()
        all_probs, all_labels = [], []
        total_loss = 0.0

        for batch in dataloader:
            batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                     for k, v in batch.items()}
            logits = self.model(
                batch["node_feats"], batch["edge_index"], batch["edge_feats"],
                batch["esm2_embedding"], batch["thermo_features"], batch["evo_features"],
                batch.get("node_mask"),
            )
            labels = batch["labels"].unsqueeze(-1)
            loss = self.criterion(logits, labels)
            total_loss += loss.item()
            all_probs.extend(torch.sigmoid(logits).squeeze(-1).cpu().tolist())
            all_labels.extend(labels.squeeze(-1).cpu().tolist())

        try:
            from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef
            import numpy as np
            labels_arr = [int(l) for l in all_labels]
            auroc = roc_auc_score(labels_arr, all_probs)
            auprc = average_precision_score(labels_arr, all_probs)
            preds = [int(p > 0.5) for p in all_probs]
            mcc = matthews_corrcoef(labels_arr, preds)
        except Exception:
            auroc = auprc = mcc = 0.0

        return {
            "loss": total_loss / max(len(all_labels), 1),
            "auroc": auroc,
            "auprc": auprc,
            "mcc": mcc,
        }

    def save_checkpoint(self, metrics: dict[str, float], tag: str = "latest") -> None:
        ckpt = {
            "model_state_dict": self.model.state_dict(),
            "ema_shadow": self.ema.shadow,
            "optimizer_state_dict": self.optimizer.state_dict(),
            "global_step": self.global_step,
            "metrics": metrics,
        }
        path = self.checkpoint_dir / f"{tag}.pt"
        torch.save(ckpt, path)
        logger.info(f"Saved checkpoint to {path} — {metrics}")

        if metrics.get("auroc", 0.0) > (1.0 - self.best_val_loss):
            self.best_val_loss = 1.0 - metrics.get("auroc", 0.0)
            best_path = self.checkpoint_dir / "best.pt"
            torch.save(ckpt, best_path)

    def train(
        self,
        train_loader: Any,
        val_loader: Any,
        n_epochs: int = 50,
        patience: int = 10,
    ) -> None:
        no_improve = 0
        for epoch in range(n_epochs):
            train_loss = self.train_epoch(train_loader, epoch)
            val_metrics = self.evaluate(val_loader)
            val_metrics["train_loss"] = train_loss

            logger.info(
                f"Epoch {epoch+1}/{n_epochs} — "
                f"train_loss={train_loss:.4f} "
                f"val_auroc={val_metrics['auroc']:.4f} "
                f"val_auprc={val_metrics['auprc']:.4f}"
            )

            self.save_checkpoint(val_metrics, tag="latest")
            current_auroc = val_metrics.get("auroc", 0.0)
            if current_auroc > (1.0 - self.best_val_loss) - 1e-4:
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= patience:
                    logger.info(f"Early stopping at epoch {epoch+1}")
                    break
