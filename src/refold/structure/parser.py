"""PDB structure parser using gemmi.

Converts AlphaFold PDB files to ProteinStructure tensors.
NaN invariant: coords[i,j,:] is NaN iff atom_mask[i,j] is False.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from refold.constants import (
    ATOM_NAMES, ATOM_TO_IDX, CA_IDX,
    AA_THREE_TO_ONE, STANDARD_AAS, N_ATOM_TYPES,
    DSSP_SS_MAP, MAX_ASA,
)
from refold.types import ProteinStructure

logger = logging.getLogger(__name__)


def parse_pdb_to_structure(
    pdb_path: Path | str,
    uniprot_id: str,
) -> ProteinStructure:
    """Parse a PDB file into a ProteinStructure tensor representation.

    Uses gemmi for PDB parsing. Only the first model and chain A (or first chain) is used.
    bfactors are treated as pLDDT values for AlphaFold structures.
    """
    import gemmi

    pdb_path = Path(pdb_path)
    structure = gemmi.read_structure(str(pdb_path))

    model = structure[0]
    chain = model[0]

    sequence_chars = []
    residue_types_list = []
    coords_list = []
    atom_mask_list = []
    bfactors_list = []

    for res in chain:
        if res.het_flag != 'A' and res.name not in AA_THREE_TO_ONE:
            continue

        res_name = res.name
        aa_one = AA_THREE_TO_ONE.get(res_name)
        if aa_one is None:
            aa_one = "X"

        from refold.constants import AA_TO_IDX
        aa_idx = AA_TO_IDX.get(aa_one, 20)

        res_coords = np.full((N_ATOM_TYPES, 3), np.nan, dtype=np.float32)
        res_mask = np.zeros(N_ATOM_TYPES, dtype=bool)
        res_bfactor = 0.0
        n_atoms = 0

        for atom in res:
            atom_name = atom.name
            if atom_name in ATOM_TO_IDX:
                idx = ATOM_TO_IDX[atom_name]
                pos = atom.pos
                res_coords[idx] = [pos.x, pos.y, pos.z]
                res_mask[idx] = True
                res_bfactor += atom.b_iso
                n_atoms += 1

        if n_atoms > 0:
            res_bfactor /= n_atoms

        sequence_chars.append(aa_one)
        residue_types_list.append(aa_idx)
        coords_list.append(res_coords)
        atom_mask_list.append(res_mask)
        bfactors_list.append(res_bfactor)

    if not sequence_chars:
        raise ValueError(f"No amino acid residues found in {pdb_path}")

    coords = np.stack(coords_list, axis=0).astype(np.float32)
    atom_mask = np.stack(atom_mask_list, axis=0)
    residue_types = np.array(residue_types_list, dtype=np.int64)
    bfactors = np.array(bfactors_list, dtype=np.float32)
    residue_mask = atom_mask[:, CA_IDX]

    # Enforce NaN invariant
    for i in range(len(sequence_chars)):
        for j in range(N_ATOM_TYPES):
            if not atom_mask[i, j]:
                coords[i, j, :] = np.nan

    sequence = "".join(sequence_chars)

    # Compute relative ASA using DSSP if available
    rel_asa: np.ndarray | None = None
    sse_ids: np.ndarray | None = None
    phi_psi: np.ndarray | None = None

    try:
        rel_asa, sse_ids, phi_psi = _compute_dssp_features(
            pdb_path, sequence, bfactors
        )
    except Exception as e:
        logger.debug(f"DSSP computation failed for {uniprot_id}: {e}")

    return ProteinStructure(
        uniprot_id=uniprot_id,
        sequence=sequence,
        coords=coords,
        residue_types=residue_types,
        residue_mask=residue_mask,
        atom_mask=atom_mask,
        bfactors=bfactors,
        sse_ids=sse_ids,
        phi_psi=phi_psi,
        rel_asa=rel_asa,
        source="alphafold",
    )


def _compute_dssp_features(
    pdb_path: Path,
    sequence: str,
    bfactors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute DSSP secondary structure, ASA, and dihedral angles."""
    try:
        from Bio.PDB import PDBParser, DSSP
        from Bio.PDB.DSSP import dssp_dict_from_pdb_file
    except ImportError:
        raise ImportError("biopython required for DSSP")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]

    try:
        dssp = DSSP(model, str(pdb_path), dssp="mkdssp")
    except Exception:
        dssp = DSSP(model, str(pdb_path))

    n_res = len(sequence)
    rel_asa = np.zeros(n_res, dtype=np.float32)
    sse_ids = np.zeros(n_res, dtype=np.int64)
    phi_psi = np.zeros((n_res, 4), dtype=np.float32)

    dssp_list = list(dssp.property_list)

    for i, record in enumerate(dssp_list):
        if i >= n_res:
            break
        ss = record[2]
        rasa = record[3]
        phi = record[4]
        psi = record[5]

        sse_ids[i] = DSSP_SS_MAP.get(ss, 7)
        rel_asa[i] = float(rasa) if rasa is not None else 0.0

        phi_rad = np.radians(phi) if phi is not None else 0.0
        psi_rad = np.radians(psi) if psi is not None else 0.0
        phi_psi[i] = [np.sin(phi_rad), np.cos(phi_rad), np.sin(psi_rad), np.cos(psi_rad)]

    return rel_asa, sse_ids, phi_psi
