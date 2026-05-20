#!/usr/bin/env python3
"""Train the rescue classifier model."""

import argparse
import logging
from pathlib import Path

import pandas as pd
import torch
from torch.utils.data import DataLoader

from refold.utils.logging import setup_logging
from refold.utils.device import get_device
from refold.constants import PROCESSED_DIR, CHECKPOINT_DIR
from refold.models.rescue_classifier.model import RescueClassifier
from refold.models.rescue_classifier.train import RescueClassifierTrainer
from refold.data.datasets.rescue_dataset import RescueDataset, collate_rescue_batch


def build_training_dataframes(
    mutations_path: Path,
    stability_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build train/val dataframes with rescue labels."""
    if not mutations_path.exists():
        raise FileNotFoundError(f"Mutations not found: {mutations_path}")

    df = pd.read_parquet(mutations_path)

    stability_df = None
    if stability_path.exists():
        stability_df = pd.read_parquet(stability_path)
        df = df.merge(
            stability_df[["uniprot_id", "position", "wildtype_aa", "mutant_aa", "ddg_kcal_mol"]],
            on=["uniprot_id", "position", "wildtype_aa", "mutant_aa"],
            how="left",
        )

    if "ddg_kcal_mol" not in df.columns:
        df["ddg_kcal_mol"] = 0.0

    df["ddg_kcal_mol"] = df["ddg_kcal_mol"].fillna(0.0)

    # Tentative rescue label: 0.5 < ΔΔG < 5.0 = positive
    df["rescue_label"] = (
        (df["ddg_kcal_mol"] > 0.5) & (df["ddg_kcal_mol"] < 5.0)
    ).astype(float)

    from refold.data.splits import create_splits
    splits = create_splits(df)

    return splits["train"], splits["val"], splits["test"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Train rescue classifier")
    parser.add_argument("--mutations", type=Path,
                        default=PROCESSED_DIR / "mutations" / "clinvar_missense.parquet")
    parser.add_argument("--stability", type=Path,
                        default=PROCESSED_DIR / "mutations" / "stability_training.parquet")
    parser.add_argument("--checkpoint-dir", type=Path,
                        default=CHECKPOINT_DIR / "rescue_classifier")
    parser.add_argument("--pretrain-epochs", type=int, default=10)
    parser.add_argument("--finetune-epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--test-run", action="store_true",
                        help="Quick test run with small dataset")
    args = parser.parse_args()

    setup_logging(args.checkpoint_dir / "logs")
    logger = logging.getLogger(__name__)

    device = torch.device(args.device) if args.device else get_device()
    logger.info(f"Training on device: {device}")

    logger.info("Building training data...")
    train_df, val_df, _ = build_training_dataframes(args.mutations, args.stability)

    if args.test_run:
        train_df = train_df.head(100)
        val_df = val_df.head(50)
        logger.info("Test run: using 100 train / 50 val samples")

    cache_dir = PROCESSED_DIR / "structures" / "rescue_cache"
    structure_dir = PROCESSED_DIR / "structures"

    train_dataset = RescueDataset(train_df, structure_dir, cache_dir / "train")
    val_dataset = RescueDataset(val_df, structure_dir, cache_dir / "val")

    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        collate_fn=collate_rescue_batch,
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        collate_fn=collate_rescue_batch,
    )

    logger.info(f"Train: {len(train_dataset):,}  Val: {len(val_dataset):,}")
    model = RescueClassifier()
    logger.info(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    trainer = RescueClassifierTrainer(
        model=model,
        device=device,
        checkpoint_dir=args.checkpoint_dir,
    )

    logger.info(f"Training for {args.finetune_epochs} epochs...")
    trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        n_epochs=args.finetune_epochs,
        patience=10,
    )
    logger.info("Training complete.")


if __name__ == "__main__":
    main()
