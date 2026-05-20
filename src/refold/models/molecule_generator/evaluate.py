"""Evaluation metrics for generated molecule libraries."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def compute_validity(smiles_list: list[Optional[str]]) -> float:
    """Fraction of non-None SMILES that parse as valid RDKit molecules."""
    try:
        from rdkit import Chem
    except ImportError:
        return 0.0

    valid = 0
    total = len(smiles_list)
    if total == 0:
        return 0.0

    for smi in smiles_list:
        if smi is None:
            continue
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                valid += 1
        except Exception:
            pass
    return valid / total


def compute_uniqueness(smiles_list: list[str]) -> float:
    """Fraction of unique canonical SMILES among valid molecules."""
    try:
        from rdkit import Chem
    except ImportError:
        return 0.0

    canonical = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical.append(Chem.MolToSmiles(mol))
        except Exception:
            pass

    if not canonical:
        return 0.0
    return len(set(canonical)) / len(canonical)


def compute_novelty(
    generated_smiles: list[str],
    training_smiles: set[str],
) -> float:
    """Fraction of generated molecules not in training set."""
    try:
        from rdkit import Chem
    except ImportError:
        return 0.0

    canonical_train = set()
    for smi in training_smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical_train.add(Chem.MolToSmiles(mol))
        except Exception:
            pass

    novel = 0
    total = 0
    for smi in generated_smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                total += 1
                if Chem.MolToSmiles(mol) not in canonical_train:
                    novel += 1
        except Exception:
            pass
    return novel / total if total > 0 else 0.0


def compute_drug_likeness_fraction(smiles_list: list[str]) -> dict[str, float]:
    """Fraction of molecules passing Lipinski RO5 and Veber filters."""
    try:
        from refold.scoring.filters import compute_lipinski_properties
        from refold.constants import MAX_MW, MAX_LOGP, MAX_HBD, MAX_HBA, MAX_ROTATABLE_BONDS, MAX_TPSA
    except ImportError:
        return {"lipinski": 0.0, "veber": 0.0}

    n_lipinski = 0
    n_veber = 0
    n_total = 0

    for smi in smiles_list:
        try:
            props = compute_lipinski_properties(smi)
            if props is None:
                continue
            n_total += 1
            if (props["mw"] <= MAX_MW and props["logp"] <= MAX_LOGP
                    and props["hbd"] <= MAX_HBD and props["hba"] <= MAX_HBA):
                n_lipinski += 1
            if props["rotatable_bonds"] <= MAX_ROTATABLE_BONDS and props["tpsa"] <= MAX_TPSA:
                n_veber += 1
        except Exception:
            pass

    if n_total == 0:
        return {"lipinski": 0.0, "veber": 0.0, "n_total": 0}

    return {
        "lipinski": n_lipinski / n_total,
        "veber": n_veber / n_total,
        "n_total": n_total,
    }


def evaluate_molecule_library(
    smiles_list: list[Optional[str]],
    training_smiles: Optional[set[str]] = None,
) -> dict[str, float]:
    """Full evaluation report for a generated molecule library."""
    valid_smiles = [s for s in smiles_list if s is not None]

    metrics: dict[str, float] = {}
    metrics["validity"] = compute_validity(smiles_list)
    metrics["uniqueness"] = compute_uniqueness(valid_smiles)

    if training_smiles is not None:
        metrics["novelty"] = compute_novelty(valid_smiles, training_smiles)

    dl = compute_drug_likeness_fraction(valid_smiles)
    metrics.update(dl)

    logger.info(f"Generation metrics: {metrics}")
    return metrics
