"""Visualization utilities for REFOLD results."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from refold.types import REFOLDResult, GeneratedMolecule

logger = logging.getLogger(__name__)


def draw_molecules_grid(
    molecules: list["GeneratedMolecule"],
    output_path: Path,
    n_cols: int = 4,
    mol_size: tuple[int, int] = (250, 200),
) -> None:
    """Draw a grid of molecules using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw

        mols = []
        labels = []
        for m in molecules:
            mol = Chem.MolFromSmiles(m.smiles)
            if mol:
                mols.append(mol)
                labels.append(
                    f"QED={m.qed_score:.2f}\nRP={m.rescue_probability:.2f}"
                )

        if not mols:
            return

        n_rows = (len(mols) + n_cols - 1) // n_cols
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=n_cols,
            subImgSize=mol_size,
            legends=labels,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        img.save(str(output_path))
        logger.info(f"Saved molecule grid to {output_path}")
    except Exception as e:
        logger.warning(f"Could not draw molecule grid: {e}")


def plot_rescue_probability_distribution(
    probabilities: list[float],
    output_path: Path,
    title: str = "Rescue Probability Distribution",
) -> None:
    """Plot histogram of rescue probabilities."""
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(probabilities, bins=20, color="steelblue", edgecolor="white")
        ax.axvline(0.5, color="red", linestyle="--", label="Threshold (0.5)")
        ax.set_xlabel("Rescue Probability")
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.legend()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        logger.warning(f"Could not plot distribution: {e}")
