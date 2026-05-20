#!/usr/bin/env python3
"""Train the pocket detector GNN model."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.constants import CHECKPOINT_DIR, PROCESSED_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="Train REFOLD pocket detector")
    parser.add_argument("--n-epochs", type=int, default=50)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--device", default=None)
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=CHECKPOINT_DIR / "pocket_detector",
    )
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)
    logger.info("Training pocket detector...")

    # Pocket detector training uses SCPDB/PDB structures with annotated pockets.
    # In absence of annotated training data, this script demonstrates the pipeline.
    logger.warning(
        "Pocket detector training requires SCPDB or similar annotated pocket data. "
        "Place pocket annotations at data/processed/pocket_annotations.parquet"
    )

    from refold.models.pocket_detector.train import train_pocket_detector

    # Stub: empty lists — replace with real data loading
    train_data: list = []
    val_data: list = []

    if not train_data:
        logger.error("No training data available. Exiting.")
        return

    metrics = train_pocket_detector(
        train_data=train_data,
        val_data=val_data,
        n_epochs=args.n_epochs,
        lr=args.lr,
        batch_size=args.batch_size,
        checkpoint_dir=args.checkpoint_dir,
        device_str=args.device,
    )
    logger.info(f"Training complete. Best epoch: {metrics.get('best_epoch')}")


if __name__ == "__main__":
    main()
