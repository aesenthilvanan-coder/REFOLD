"""Chemistry utility functions."""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


def xyz_to_smiles(atom_positions: "np.ndarray", atom_types: list[str]) -> str | None:
    """Convert 3D atom coordinates to SMILES via Open Babel."""
    import numpy as np

    xyz_lines = [f"{len(atom_types)}", "REFOLD generated molecule"]
    for atom, pos in zip(atom_types, atom_positions):
        xyz_lines.append(f"{atom}  {pos[0]:.4f}  {pos[1]:.4f}  {pos[2]:.4f}")
    xyz_content = "\n".join(xyz_lines)

    with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w", delete=False) as f:
        f.write(xyz_content)
        xyz_path = f.name

    try:
        result = subprocess.run(
            ["obabel", xyz_path, "-osmi", "--gen3D"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            smiles = result.stdout.strip().split()[0]
            return smiles if smiles else None
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.debug(f"obabel failed: {e}")
        return None
    finally:
        Path(xyz_path).unlink(missing_ok=True)


def sanitize_smiles(smiles: str) -> str | None:
    """Sanitize and canonicalize a SMILES string using RDKit."""
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def kcal_to_um(delta_g_kcal: float, temperature_k: float = 298.15) -> float:
    """Convert binding free energy (kcal/mol) to Kd in μM."""
    import math

    R = 1.987e-3  # kcal/mol/K
    kd_molar = math.exp(delta_g_kcal / (R * temperature_k))
    return kd_molar * 1e6
