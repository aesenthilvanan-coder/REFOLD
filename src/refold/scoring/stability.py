"""Binding affinity prediction and stabilization scoring."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from refold.types import Pocket, GeneratedMolecule

logger = logging.getLogger(__name__)


def predict_binding_affinity_vina(
    smiles: str,
    pocket: "Pocket",
    protein_pdb_path: Path,
) -> float | None:
    """Run AutoDock Vina docking. Returns ΔG in kcal/mol or None if unavailable."""
    if shutil.which("vina") is None:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)

            ligand_path = Path(tmp) / "ligand.mol2"
            Chem.MolToMolFile(mol, str(ligand_path))

            center = pocket.center
            config_path = Path(tmp) / "vina.conf"
            with open(config_path, "w") as f:
                f.write(f"receptor = {protein_pdb_path}\n")
                f.write(f"ligand = {ligand_path}\n")
                f.write(f"center_x = {center[0]:.2f}\n")
                f.write(f"center_y = {center[1]:.2f}\n")
                f.write(f"center_z = {center[2]:.2f}\n")
                f.write("size_x = 20\nsize_y = 20\nsize_z = 20\n")
                f.write("out = result.pdbqt\nexhaust = 8\nnum_modes = 1\n")

            result = subprocess.run(
                ["vina", "--config", str(config_path)],
                capture_output=True, text=True, timeout=120, cwd=tmp,
            )

            for line in result.stdout.split("\n"):
                if line.strip().startswith("1"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return float(parts[1])
        except Exception as e:
            logger.debug(f"Vina docking failed: {e}")

    return None


def predict_binding_affinity_ml(
    molecule: "GeneratedMolecule",
    pocket: "Pocket",
) -> float:
    """Fast heuristic binding affinity estimate in kcal/mol."""
    druggability_bonus = pocket.druggability_score * 3.0
    hydrophobic_bonus = pocket.hydrophobicity * 1.5

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(molecule.smiles)
        mw = Descriptors.ExactMolWt(mol) if mol else 300.0
    except Exception:
        mw = 300.0

    size_penalty = max(0.0, (mw - 300.0) / 200.0)

    affinity_kcal = -(druggability_bonus + hydrophobic_bonus - size_penalty)
    affinity_kcal = max(-15.0, min(-1.0, affinity_kcal))
    return float(affinity_kcal)


def compute_stabilization_score(
    affinity_kcal: float,
    pocket: "Pocket",
    ddg_protein: float,
) -> float:
    """Compute rescue stabilization score.

    ΔG_rescue = |ΔG_bind| × detection_freq × min(ΔΔG/5, 1)
    """
    if affinity_kcal >= 0:
        return 0.0
    ddg_factor = min(abs(ddg_protein) / 5.0, 1.0)
    score = abs(affinity_kcal) * pocket.detection_frequency * ddg_factor
    return float(score)
