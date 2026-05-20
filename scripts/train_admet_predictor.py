#!/usr/bin/env python3
"""Train the multi-task ADMET predictor."""

import argparse
import logging
from pathlib import Path

import numpy as np

from refold.utils.logging import setup_logging
from refold.constants import CHECKPOINT_DIR, PROCESSED_DIR


def load_admet_training_data() -> tuple:
    """Load ADMET training data from processed ChEMBL/BindingDB files."""
    chembl_path = PROCESSED_DIR / "chembl_binding.parquet"
    bindingdb_path = PROCESSED_DIR / "bindingdb_binding.parquet"

    import pandas as pd
    from refold.models.admet_predictor.evaluate import ADMET_TASK_NAMES

    if not chembl_path.exists() and not bindingdb_path.exists():
        return [], np.empty((0, len(ADMET_TASK_NAMES))), [], np.empty((0, len(ADMET_TASK_NAMES)))

    dfs = []
    if chembl_path.exists():
        dfs.append(pd.read_parquet(chembl_path))
    if bindingdb_path.exists():
        dfs.append(pd.read_parquet(bindingdb_path))

    df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    if df.empty or "smiles" not in df.columns:
        return [], np.empty((0, len(ADMET_TASK_NAMES))), [], np.empty((0, len(ADMET_TASK_NAMES)))

    smiles = df["smiles"].dropna().tolist()
    # Stub labels — real implementation would source task-specific labels
    n = len(smiles)
    n_tasks = len(ADMET_TASK_NAMES)
    labels = np.full((n, n_tasks), np.nan, dtype=np.float32)

    # 70/30 train/val split
    split = int(0.7 * n)
    return (
        smiles[:split], labels[:split],
        smiles[split:], labels[split:],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Train REFOLD ADMET predictor")
    parser.add_argument("--n-epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--device", default=None)
    parser.add_argument(
        "--checkpoint-dir", type=Path,
        default=CHECKPOINT_DIR / "admet_predictor",
    )
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)
    logger.info("Loading ADMET training data...")

    train_smiles, train_labels, val_smiles, val_labels = load_admet_training_data()

    if not train_smiles:
        logger.error(
            "No ADMET training data. Run 'make download' and 'make preprocess' first."
        )
        return

    logger.info(f"Training on {len(train_smiles)} molecules...")

    from refold.models.admet_predictor.train import train_admet_predictor
    metrics = train_admet_predictor(
        train_smiles=train_smiles,
        train_labels=train_labels,
        val_smiles=val_smiles,
        val_labels=val_labels,
        n_epochs=args.n_epochs,
        lr=args.lr,
        batch_size=args.batch_size,
        checkpoint_dir=args.checkpoint_dir,
        device_str=args.device,
    )
    logger.info(f"Training complete. Best epoch: {metrics.get('best_epoch')}")


if __name__ == "__main__":
    main()
