"""Core data types for REFOLD."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np


class MutationClass(Enum):
    MISFOLDING = "MISFOLDING"
    FUNCTIONAL_DISRUPTION = "FUNCTIONAL_DISRUPTION"
    TRAFFICKING = "TRAFFICKING"
    AGGREGATION = "AGGREGATION"
    UNKNOWN = "UNKNOWN"


class RescueAmenability(Enum):
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    UNRESCUABLE = "UNRESCUABLE"


class PocketType(Enum):
    WILDTYPE_ORTHOSTERIC = "WILDTYPE_ORTHOSTERIC"
    WILDTYPE_ALLOSTERIC = "WILDTYPE_ALLOSTERIC"
    TRANSIENT_MISFOLDING = "TRANSIENT_MISFOLDING"
    CRYPTIC = "CRYPTIC"


@dataclass
class Mutation:
    uniprot_id: str
    position: int          # 1-based residue index
    wildtype_aa: str       # single-letter code
    mutant_aa: str         # single-letter code
    gene_name: str = ""
    disease: str = ""
    clinvar_id: str = ""
    source: str = ""

    @property
    def hgvs(self) -> str:
        """Return HGVS protein notation p.{WT}{pos}{MT}."""
        return f"p.{self.wildtype_aa}{self.position}{self.mutant_aa}"

    def __str__(self) -> str:
        return f"{self.gene_name or self.uniprot_id}:{self.hgvs}"

    def to_key(self) -> str:
        """Canonical string key for deduplication and benchmark lookup."""
        return f"{self.uniprot_id}_{self.position}_{self.wildtype_aa}_{self.mutant_aa}"


@dataclass
class ProteinStructure:
    uniprot_id: str
    sequence: str                        # single-letter AA sequence, length N_res
    coords: np.ndarray                   # [N_res, 37, 3] float32 — NaN where atom absent
    residue_types: np.ndarray            # [N_res] int64, 0-19=standard AA, 20=UNK
    residue_mask: np.ndarray             # [N_res] bool — True for resolved residues
    atom_mask: np.ndarray                # [N_res, 37] bool
    bfactors: np.ndarray                 # [N_res] float32 (pLDDT from AF, 0-100)
    sse_ids: np.ndarray | None = None    # [N_res] int8 (DSSP_SS_MAP encoding)
    phi_psi: np.ndarray | None = None   # [N_res, 4] float32: sin_phi, cos_phi, sin_psi, cos_psi
    rel_asa: np.ndarray | None = None   # [N_res] float32, relative solvent accessibility
    source: str = "alphafold"

    def __post_init__(self) -> None:
        n = len(self.sequence)
        assert self.coords.shape == (n, 37, 3), (
            f"coords shape {self.coords.shape} != ({n}, 37, 3)"
        )
        assert self.residue_types.shape == (n,), (
            f"residue_types shape {self.residue_types.shape} != ({n},)"
        )
        assert self.atom_mask.shape == (n, 37), (
            f"atom_mask shape {self.atom_mask.shape} != ({n}, 37)"
        )
        # NaN invariant: coords[i,j,:] is NaN iff atom_mask[i,j] is False
        masked_not_nan = (~self.atom_mask) & (~np.any(np.isnan(self.coords), axis=-1))
        assert not np.any(masked_not_nan), "atom_mask/coords NaN invariant violated"

    @property
    def n_residues(self) -> int:
        return len(self.sequence)

    @property
    def ca_coords(self) -> np.ndarray:
        """[N_res, 3] Cα coordinates; NaN for unresolved residues."""
        return self.coords[:, 1, :]  # CA_IDX = 1

    @property
    def plddt(self) -> np.ndarray:
        """Alias for bfactors (pLDDT values from AlphaFold)."""
        return self.bfactors

    @property
    def disordered_mask(self) -> np.ndarray:
        """[N_res] bool — True where pLDDT < 50 (highly disordered)."""
        return self.bfactors < 50.0


@dataclass
class MutantStructure:
    wildtype: ProteinStructure
    mutation: Mutation
    mutant_coords: np.ndarray            # [N_res, 37, 3]
    mutant_residue_types: np.ndarray     # [N_res] int64
    ddg_predicted: float                 # kcal/mol, positive = destabilising
    ddg_esm1v: float | None = None
    ddg_gnn: float | None = None
    relaxed: bool = False


@dataclass
class Pocket:
    pocket_id: str
    pocket_type: PocketType
    center: np.ndarray            # [3] float32, Å
    volume: float                 # Å³
    druggability_score: float     # [0, 1]
    hydrophobicity: float         # [0, 1]
    residue_indices: list[int]    # 0-based indices into protein sequence
    alpha_sphere_coords: np.ndarray  # [n_spheres, 3] float32
    detection_frequency: float = 1.0   # fraction of ensemble conformations where detected

    @property
    def is_transient(self) -> bool:
        return self.pocket_type == PocketType.TRANSIENT_MISFOLDING

    @property
    def is_druggable(self) -> bool:
        from refold.constants import DRUGGABILITY_THRESHOLD
        return self.druggability_score >= DRUGGABILITY_THRESHOLD

    def __str__(self) -> str:
        return (
            f"Pocket({self.pocket_id}, type={self.pocket_type.value}, "
            f"vol={self.volume:.0f}Å³, drug={self.druggability_score:.2f}, "
            f"freq={self.detection_frequency:.2f})"
        )


@dataclass
class GeneratedMolecule:
    smiles: str
    pocket_id: str
    predicted_affinity_kcal: float = 0.0
    sa_score: float = 0.0
    qed_score: float = 0.0
    rescue_probability: float = 0.0

    # Lipinski properties
    mw: float = 0.0
    logp: float = 0.0
    hbd: int = 0
    hba: int = 0
    rotatable_bonds: int = 0
    tpsa: float = 0.0
    passes_lipinski: bool = False
    passes_veber: bool = False

    # Additional filters
    is_pains: bool = False
    passes_all_filters: bool = False

    # ADMET predictions
    herg_inhibition_prob: float = 0.0
    bbb_permeability: float = 0.0
    cyp3a4_inhibition_prob: float = 0.0
    oral_bioavailability: float = 0.0
    solubility_logS: float = 0.0

    # Metadata
    generation_step: int = 0
    rank: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "smiles": self.smiles,
            "pocket_id": self.pocket_id,
            "predicted_affinity_kcal": self.predicted_affinity_kcal,
            "sa_score": self.sa_score,
            "qed_score": self.qed_score,
            "rescue_probability": self.rescue_probability,
            "mw": self.mw,
            "logp": self.logp,
            "hbd": self.hbd,
            "hba": self.hba,
            "rotatable_bonds": self.rotatable_bonds,
            "tpsa": self.tpsa,
            "passes_lipinski": self.passes_lipinski,
            "passes_veber": self.passes_veber,
            "is_pains": self.is_pains,
            "passes_all_filters": self.passes_all_filters,
            "herg_inhibition_prob": self.herg_inhibition_prob,
            "bbb_permeability": self.bbb_permeability,
            "cyp3a4_inhibition_prob": self.cyp3a4_inhibition_prob,
            "oral_bioavailability": self.oral_bioavailability,
            "solubility_logS": self.solubility_logS,
            "rank": self.rank,
        }


@dataclass
class REFOLDResult:
    mutation: Mutation
    mutation_class: MutationClass
    rescue_amenability: RescueAmenability
    rescue_amenability_prob: float
    rescue_probability: float
    ddg_predicted: float
    n_pockets_detected: int
    n_molecules_generated: int
    n_molecules_passing_filters: int
    top_candidates: list[GeneratedMolecule] = field(default_factory=list)
    pockets: list[Pocket] = field(default_factory=list)
    all_molecules: list[GeneratedMolecule] = field(default_factory=list)
    error_message: str | None = None
    runtime_seconds: float = 0.0

    @property
    def is_rescue_amenable(self) -> bool:
        return self.rescue_amenability in (RescueAmenability.HIGH, RescueAmenability.MODERATE)

    @property
    def has_druggable_pocket(self) -> bool:
        return any(p.is_druggable for p in self.pockets)

    def to_dict(self) -> dict[str, Any]:
        return {
            "mutation": {
                "uniprot_id": self.mutation.uniprot_id,
                "gene_name": self.mutation.gene_name,
                "position": self.mutation.position,
                "wildtype_aa": self.mutation.wildtype_aa,
                "mutant_aa": self.mutation.mutant_aa,
                "hgvs": self.mutation.hgvs,
                "disease": self.mutation.disease,
                "clinvar_id": self.mutation.clinvar_id,
            },
            "classification": {
                "mutation_class": self.mutation_class.value,
                "rescue_amenability": self.rescue_amenability.value,
                "rescue_amenability_prob": self.rescue_amenability_prob,
                "rescue_probability": self.rescue_probability,
                "ddg_predicted": self.ddg_predicted,
            },
            "pockets": [
                {
                    "pocket_id": p.pocket_id,
                    "pocket_type": p.pocket_type.value,
                    "volume": p.volume,
                    "druggability_score": p.druggability_score,
                    "detection_frequency": p.detection_frequency,
                    "n_residues": len(p.residue_indices),
                }
                for p in self.pockets
            ],
            "molecules": {
                "n_generated": self.n_molecules_generated,
                "n_passing_filters": self.n_molecules_passing_filters,
                "top_candidates": [m.to_dict() for m in self.top_candidates],
            },
            "runtime_seconds": self.runtime_seconds,
            "error": self.error_message,
        }
