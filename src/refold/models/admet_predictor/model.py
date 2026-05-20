"""ADMET property predictor using molecular fingerprints."""

from __future__ import annotations

import logging

import torch
import torch.nn as nn

from refold.constants import ADMET_HIDDEN_DIMS, ADMET_N_TASKS

logger = logging.getLogger(__name__)

ADMET_TASK_NAMES = [
    "bbb_permeability",
    "hia_absorption",
    "pgp_substrate",
    "herg_inhibition",
    "ames_toxicity",
    "ld50_log",
    "cyp3a4_inhibition",
    "solubility_logS",
]


class ADMETPredictor(nn.Module):
    """Multi-task ADMET property predictor.

    Takes Morgan fingerprints and predicts 8 ADMET properties.
    Tasks: BBB, HIA, Pgp, hERG, AMES, LD50, CYP3A4, solubility.
    """

    def __init__(
        self,
        input_dim: int = 2048,
        hidden_dims: list[int] = ADMET_HIDDEN_DIMS,
        n_tasks: int = ADMET_N_TASKS,
        dropout: float = 0.2,
    ):
        super().__init__()

        layers: list[nn.Module] = []
        in_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.GELU(),
                nn.Dropout(dropout),
            ])
            in_dim = hidden_dim

        self.encoder = nn.Sequential(*layers)
        self.task_heads = nn.ModuleList([
            nn.Linear(in_dim, 1) for _ in range(n_tasks)
        ])

    def forward(self, fingerprints: torch.Tensor) -> torch.Tensor:
        """Forward pass.

        Args:
            fingerprints: [B, input_dim] Morgan fingerprints.

        Returns:
            [B, n_tasks] predictions (raw logits for binary tasks).
        """
        encoded = self.encoder(fingerprints)
        outputs = [head(encoded) for head in self.task_heads]
        return torch.cat(outputs, dim=-1)

    def predict_properties(self, smiles: str) -> dict[str, float]:
        """Predict ADMET properties for a SMILES string."""
        try:
            fp = _smiles_to_fingerprint(smiles)
            if fp is None:
                return {task: 0.5 for task in ADMET_TASK_NAMES}

            device = next(self.parameters()).device
            fp_tensor = torch.tensor(fp, dtype=torch.float32).unsqueeze(0).to(device)

            with torch.no_grad():
                logits = self(fp_tensor).squeeze(0)
                probs = torch.sigmoid(logits).cpu().tolist()

            return {name: float(p) for name, p in zip(ADMET_TASK_NAMES, probs)}
        except Exception as e:
            logger.debug(f"ADMET prediction failed for {smiles}: {e}")
            return {task: 0.5 for task in ADMET_TASK_NAMES}


def _smiles_to_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> list[int] | None:
    """Convert SMILES to Morgan fingerprint."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return list(fp)
    except Exception:
        return None
