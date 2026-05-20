"""PyTorch dataset for molecule generator training."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from refold.constants import (
    MOL_ATOM_TYPES, N_MOL_ATOM_TYPES, N_ATOMS_MIN, N_ATOMS_MAX,
    POCKET_NODE_FEAT_DIM, POCKET_SCALAR_DIM,
)

logger = logging.getLogger(__name__)


class MoleculeDataset:
    """
    Dataset of (pocket, molecule) pairs for diffusion model training.

    Each item:
      - pocket_node_feats: [K, POCKET_NODE_FEAT_DIM] float32
      - pocket_scalars: [1, POCKET_SCALAR_DIM] float32
      - atom_positions: [N_atoms, 3] float32
      - atom_types: [N_atoms] int64
      - n_atoms: int
    """

    def __init__(
        self,
        pockets: list,
        molecules_by_pocket: dict[str, list],
        max_atoms: int = N_ATOMS_MAX,
    ):
        self.pockets = pockets
        self.molecules_by_pocket = molecules_by_pocket
        self.max_atoms = max_atoms
        self._items = self._build_index()

    def _build_index(self) -> list[tuple[int, int]]:
        """Build flat (pocket_idx, mol_idx) index."""
        items = []
        for pi, pocket in enumerate(self.pockets):
            mols = self.molecules_by_pocket.get(pocket.pocket_id, [])
            for mi in range(len(mols)):
                items.append((pi, mi))
        return items

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, idx: int) -> dict:
        pi, mi = self._items[idx]
        pocket = self.pockets[pi]
        mols = self.molecules_by_pocket.get(pocket.pocket_id, [])

        if mi >= len(mols):
            return self._dummy_item()

        try:
            return self._build_item(pocket, mols[mi])
        except Exception as e:
            logger.debug(f"MoleculeDataset item {idx} failed: {e}")
            return self._dummy_item()

    def _build_item(self, pocket, molecule) -> dict:
        from refold.models.molecule_generator.encoder import (
            build_pocket_node_features, build_pocket_scalar_features,
        )

        pocket_node_feats = build_pocket_node_features(pocket)
        pocket_scalars = build_pocket_scalar_features(pocket)

        # Molecule 3D representation
        if hasattr(molecule, "atom_positions") and molecule.atom_positions is not None:
            atom_positions = np.array(molecule.atom_positions, dtype=np.float32)
            atom_types = np.array(molecule.atom_types, dtype=np.int64)
        else:
            # Fallback: generate 3D coords from SMILES
            atom_positions, atom_types = _smiles_to_3d(molecule.smiles)

        n_atoms = len(atom_positions)
        # Pad / truncate to max_atoms
        if n_atoms > self.max_atoms:
            atom_positions = atom_positions[:self.max_atoms]
            atom_types = atom_types[:self.max_atoms]
            n_atoms = self.max_atoms

        return {
            "pocket_node_feats": pocket_node_feats,
            "pocket_scalars": pocket_scalars,
            "atom_positions": atom_positions,
            "atom_types": atom_types,
            "n_atoms": n_atoms,
        }

    def _dummy_item(self) -> dict:
        n_sphere = 5
        return {
            "pocket_node_feats": np.zeros((n_sphere, POCKET_NODE_FEAT_DIM), dtype=np.float32),
            "pocket_scalars": np.zeros((1, POCKET_SCALAR_DIM), dtype=np.float32),
            "atom_positions": np.zeros((N_ATOMS_MIN, 3), dtype=np.float32),
            "atom_types": np.zeros(N_ATOMS_MIN, dtype=np.int64),
            "n_atoms": N_ATOMS_MIN,
        }


def _smiles_to_3d(smiles: str) -> tuple[np.ndarray, np.ndarray]:
    """Generate 3D coords from SMILES using RDKit ETKDG."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)
        conf = mol.GetConformer()
        positions = np.array(conf.GetPositions(), dtype=np.float32)
        types = []
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in MOL_ATOM_TYPES:
                types.append(MOL_ATOM_TYPES.index(sym))
            else:
                types.append(0)  # default to Carbon
        return positions, np.array(types, dtype=np.int64)
    except Exception:
        n = N_ATOMS_MIN
        return np.random.randn(n, 3).astype(np.float32), np.zeros(n, dtype=np.int64)
