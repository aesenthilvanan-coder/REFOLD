#!/usr/bin/env python3
"""Train the pocket-conditioned molecule generator."""

import argparse
import logging
from pathlib import Path

import torch
from torch.utils.data import DataLoader, TensorDataset

from refold.utils.logging import setup_logging
from refold.utils.device import get_device
from refold.constants import CHECKPOINT_DIR, N_MOL_ATOM_TYPES, N_ATOMS_MAX
from refold.models.molecule_generator.diffusion import REFOLDDiffusionModel
from refold.models.molecule_generator.train import MoleculeGeneratorTrainer


def make_synthetic_loader(n_samples: int, batch_size: int) -> DataLoader:
    """Create a synthetic dataloader for testing without real data."""
    import numpy as np
    from refold.constants import POCKET_NODE_FEAT_DIM, POCKET_SCALAR_DIM

    n_atoms = 15
    n_spheres = 25

    data = {
        "coords": torch.randn(n_samples, n_atoms, 3),
        "atom_types": torch.zeros(n_samples, n_atoms, N_MOL_ATOM_TYPES),
        "pocket_sphere_feats": torch.randn(n_samples, n_spheres, POCKET_NODE_FEAT_DIM),
        "pocket_scalars": torch.randn(n_samples, POCKET_SCALAR_DIM),
    }
    data["atom_types"][:, :, 0] = 1.0

    dataset = list(
        {k: v[i] for k, v in data.items()}
        for i in range(n_samples)
    )

    def collate(batch):
        return {k: torch.stack([item[k] for item in batch]) for k in batch[0]}

    return DataLoader(dataset, batch_size=batch_size, shuffle=True, collate_fn=collate)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train molecule generator")
    parser.add_argument("--checkpoint-dir", type=Path,
                        default=CHECKPOINT_DIR / "molecule_generator")
    parser.add_argument("--n-steps", type=int, default=100_000)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--T", type=int, default=1000)
    parser.add_argument("--device", default=None)
    parser.add_argument("--synthetic-data", action="store_true",
                        help="Use synthetic data for testing")
    args = parser.parse_args()

    setup_logging(args.checkpoint_dir / "logs")
    logger = logging.getLogger(__name__)

    device = torch.device(args.device) if args.device else get_device()
    logger.info(f"Training on device: {device}")

    if args.synthetic_data:
        logger.info("Using synthetic training data")
        train_loader = make_synthetic_loader(1000, args.batch_size)
        val_loader = make_synthetic_loader(100, args.batch_size)
    else:
        raise NotImplementedError(
            "Real molecule dataset not yet implemented. "
            "Use --synthetic-data for testing, or provide BindingDB pocket/molecule pairs."
        )

    model = REFOLDDiffusionModel(T=args.T)
    logger.info(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    trainer = MoleculeGeneratorTrainer(
        model=model,
        device=device,
        checkpoint_dir=args.checkpoint_dir,
        batch_size=args.batch_size,
    )

    trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        n_steps=args.n_steps,
    )
    logger.info("Training complete.")


if __name__ == "__main__":
    main()
