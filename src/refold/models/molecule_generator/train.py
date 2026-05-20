"""Molecule generator training loop."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn

from refold.constants import (
    CHECKPOINT_DIR, EMA_DECAY, GRAD_CLIP_NORM,
)
from refold.models.molecule_generator.diffusion import REFOLDDiffusionModel
from refold.models.rescue_classifier.train import EMA

logger = logging.getLogger(__name__)


class MoleculeGeneratorTrainer:
    """Training loop for the REFOLD diffusion model."""

    def __init__(
        self,
        model: REFOLDDiffusionModel,
        device: torch.device,
        checkpoint_dir: Path | None = None,
        learning_rate: float = 1e-4,
        weight_decay: float = 0.01,
        gradient_accumulation: int = 8,
    ):
        self.model = model.to(device)
        self.device = device
        self.checkpoint_dir = checkpoint_dir or CHECKPOINT_DIR / "molecule_generator"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=learning_rate, weight_decay=weight_decay
        )
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=100_000, eta_min=1e-6
        )
        self.ema = EMA(model, decay=EMA_DECAY)
        self.gradient_accumulation = gradient_accumulation
        self.best_val_loss = float("inf")
        self.global_step = 0

    def train_step(self, batch: dict[str, torch.Tensor]) -> float:
        """Single training step with gradient accumulation."""
        batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}

        B = batch["coords"].shape[0]
        t = torch.randint(0, self.model.T, (B,), device=self.device)

        loss = self.model(
            batch["coords"],
            batch["atom_types"],
            batch["pocket_sphere_feats"],
            batch["pocket_scalars"],
            t,
        )
        loss = loss / self.gradient_accumulation
        loss.backward()

        scaled_loss = loss.item() * self.gradient_accumulation

        if (self.global_step + 1) % self.gradient_accumulation == 0:
            nn.utils.clip_grad_norm_(self.model.parameters(), GRAD_CLIP_NORM)
            self.optimizer.step()
            self.scheduler.step()
            self.optimizer.zero_grad()
            self.ema.update(self.model)

            # MPS memory cleanup every 200 steps
            if self.global_step % 200 == 0 and self.device.type == "mps":
                import torch
                torch.mps.empty_cache()

        self.global_step += 1
        return scaled_loss

    @torch.no_grad()
    def _validate(self, val_loader: Any, n_batches: int = 50) -> float:
        self.model.eval()
        total_loss = 0.0
        count = 0
        for i, batch in enumerate(val_loader):
            if i >= n_batches:
                break
            batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                     for k, v in batch.items()}
            B = batch["coords"].shape[0]
            t = torch.randint(0, self.model.T, (B,), device=self.device)
            loss = self.model(
                batch["coords"], batch["atom_types"],
                batch["pocket_sphere_feats"], batch["pocket_scalars"], t,
            )
            total_loss += loss.item()
            count += 1
        self.model.train()
        return total_loss / max(count, 1)

    def _save_checkpoint(self, val_loss: float, tag: str = "latest") -> None:
        ckpt = {
            "model_state_dict": self.model.state_dict(),
            "ema_shadow": self.ema.shadow,
            "optimizer_state_dict": self.optimizer.state_dict(),
            "scheduler_state_dict": self.scheduler.state_dict(),
            "global_step": self.global_step,
            "val_loss": val_loss,
        }
        torch.save(ckpt, self.checkpoint_dir / f"{tag}.pt")
        if val_loss < self.best_val_loss:
            self.best_val_loss = val_loss
            torch.save(ckpt, self.checkpoint_dir / "best.pt")
            logger.info(f"New best model at step {self.global_step}, val_loss={val_loss:.4f}")

    def train(
        self,
        train_loader: Any,
        val_loader: Any | None = None,
        n_steps: int = 100_000,
        val_every: int = 1000,
        log_every: int = 100,
    ) -> None:
        self.model.train()
        self.optimizer.zero_grad()

        running_loss = 0.0
        data_iter = iter(train_loader)

        for step in range(n_steps):
            try:
                batch = next(data_iter)
            except StopIteration:
                data_iter = iter(train_loader)
                batch = next(data_iter)

            loss = self.train_step(batch)
            running_loss += loss

            if (step + 1) % log_every == 0:
                avg_loss = running_loss / log_every
                lr = self.optimizer.param_groups[0]["lr"]
                logger.info(
                    f"Step {step+1}/{n_steps} — loss={avg_loss:.4f} lr={lr:.2e}"
                )
                running_loss = 0.0

            if val_loader and (step + 1) % val_every == 0:
                val_loss = self._validate(val_loader)
                self._save_checkpoint(val_loss, tag="latest")
                logger.info(f"Validation at step {step+1}: val_loss={val_loss:.4f}")

        self._save_checkpoint(self.best_val_loss, tag="latest")
        logger.info(f"Training complete after {n_steps} steps")
