"""Molecular fingerprint features for ADMET prediction."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

_FP_SIZE = 2048
_MORGAN_RADIUS = 2


def smiles_to_morgan_fp(smiles: str, radius: int = _MORGAN_RADIUS, n_bits: int = _FP_SIZE) -> Optional[np.ndarray]:
    """Convert SMILES to Morgan fingerprint bit vector."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp, dtype=np.float32)
    except Exception as e:
        logger.debug(f"Morgan FP failed for {smiles}: {e}")
        return None


def smiles_to_rdkit_fp(smiles: str, n_bits: int = 2048) -> Optional[np.ndarray]:
    """Convert SMILES to RDKit topological fingerprint."""
    try:
        from rdkit import Chem
        from rdkit.Chem import RDKFingerprint
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = RDKFingerprint(mol, fpSize=n_bits)
        return np.array(fp, dtype=np.float32)
    except Exception as e:
        logger.debug(f"RDKit FP failed for {smiles}: {e}")
        return None


def smiles_to_maccs_keys(smiles: str) -> Optional[np.ndarray]:
    """Convert SMILES to MACCS keys (167 bits)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import MACCSkeys
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        keys = MACCSkeys.GenMACCSKeys(mol)
        return np.array(keys, dtype=np.float32)
    except Exception as e:
        logger.debug(f"MACCS keys failed for {smiles}: {e}")
        return None


def build_admet_features(smiles: str) -> Optional[np.ndarray]:
    """
    Build combined molecular feature vector for ADMET prediction.

    Concatenates Morgan FP (2048) + MACCS (167) = 2215-dim float32.
    Returns None if RDKit cannot parse the SMILES.
    """
    morgan = smiles_to_morgan_fp(smiles)
    maccs = smiles_to_maccs_keys(smiles)

    if morgan is None or maccs is None:
        return None

    return np.concatenate([morgan, maccs], axis=0).astype(np.float32)


ADMET_FEAT_DIM: int = _FP_SIZE + 167  # 2215
