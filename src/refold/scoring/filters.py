"""Drug-likeness filters: Lipinski RO5, Veber, SA score, QED, PAINS."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from refold.types import GeneratedMolecule

from refold.constants import (
    MAX_MW, MAX_LOGP, MAX_HBD, MAX_HBA,
    MAX_ROTATABLE_BONDS, MAX_TPSA, MAX_SA_SCORE, MIN_QED,
)

logger = logging.getLogger(__name__)


def compute_lipinski_properties(smiles: str) -> dict:
    """Compute Lipinski RO5 and Veber properties from SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return _empty_lipinski()

        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        passes_lipinski = (
            mw <= MAX_MW and logp <= MAX_LOGP and hbd <= MAX_HBD and hba <= MAX_HBA
        )
        passes_veber = rotatable <= MAX_ROTATABLE_BONDS and tpsa <= MAX_TPSA

        return {
            "mw": mw, "logp": logp, "hbd": hbd, "hba": hba,
            "rotatable_bonds": rotatable, "tpsa": tpsa,
            "passes_lipinski": passes_lipinski,
            "passes_veber": passes_veber,
        }
    except Exception as e:
        logger.debug(f"Lipinski computation failed for {smiles}: {e}")
        return _empty_lipinski()


def _empty_lipinski() -> dict:
    return {
        "mw": 0.0, "logp": 0.0, "hbd": 0, "hba": 0,
        "rotatable_bonds": 0, "tpsa": 0.0,
        "passes_lipinski": False, "passes_veber": False,
    }


def compute_sa_score(smiles: str) -> float:
    """Compute synthetic accessibility score (1=easy, 10=hard)."""
    try:
        from rdkit import Chem
        try:
            from rdkit.Contrib.SA_Score import sascorer
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 10.0
            return sascorer.calculateScore(mol)
        except ImportError:
            pass

        # Fallback heuristic based on ring count and stereocenters
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 10.0
        from rdkit.Chem import rdMolDescriptors
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        try:
            n_stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        except AttributeError:
            n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        score = 1.5 + 0.4 * n_rings + 0.7 * n_stereo + 1.0 * n_spiro
        return float(min(score, 10.0))
    except Exception:
        return 5.0


def compute_qed(smiles: str) -> float:
    """Compute QED (Quantitative Estimate of Drug-likeness) score [0,1]."""
    try:
        from rdkit import Chem
        from rdkit.Chem import QED

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
        return float(QED.qed(mol))
    except Exception:
        return 0.0


def check_pains(smiles: str) -> bool:
    """Check if molecule matches PAINS filters. Returns True if PAINS positive."""
    try:
        from rdkit import Chem
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return True

        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        return catalog.HasMatch(mol)
    except Exception:
        return False


def apply_all_filters(molecule: "GeneratedMolecule") -> "GeneratedMolecule":
    """Populate all scoring fields on a GeneratedMolecule in-place."""
    lip = compute_lipinski_properties(molecule.smiles)
    molecule.mw = lip["mw"]
    molecule.logp = lip["logp"]
    molecule.hbd = lip["hbd"]
    molecule.hba = lip["hba"]
    molecule.rotatable_bonds = lip["rotatable_bonds"]
    molecule.tpsa = lip["tpsa"]
    molecule.passes_lipinski = lip["passes_lipinski"]
    molecule.passes_veber = lip["passes_veber"]

    molecule.sa_score = compute_sa_score(molecule.smiles)
    molecule.qed_score = compute_qed(molecule.smiles)
    molecule.is_pains = check_pains(molecule.smiles)

    molecule.passes_all_filters = (
        molecule.passes_lipinski
        and molecule.passes_veber
        and molecule.sa_score <= MAX_SA_SCORE
        and molecule.qed_score >= MIN_QED
        and not molecule.is_pains
    )
    return molecule


def filter_molecule_library(
    molecules: list["GeneratedMolecule"],
    require_lipinski: bool = True,
    require_veber: bool = True,
    max_sa_score: float = MAX_SA_SCORE,
    min_qed: float = MIN_QED,
    remove_pains: bool = True,
) -> list["GeneratedMolecule"]:
    """Filter and sort molecule library. Returns molecules sorted by QED descending."""
    filtered = []
    for mol in molecules:
        apply_all_filters(mol)

        if require_lipinski and not mol.passes_lipinski:
            continue
        if require_veber and not mol.passes_veber:
            continue
        if mol.sa_score > max_sa_score:
            continue
        if mol.qed_score < min_qed:
            continue
        if remove_pains and mol.is_pains:
            continue
        filtered.append(mol)

    filtered.sort(key=lambda m: m.qed_score, reverse=True)
    for i, mol in enumerate(filtered):
        mol.rank = i + 1
    return filtered
