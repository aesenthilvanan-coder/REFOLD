"""DSSP-based secondary structure and dihedral angle features."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from refold.constants import DSSP_SS_MAP, CA_IDX, N_IDX, C_IDX

logger = logging.getLogger(__name__)


def compute_backbone_dihedrals(coords: np.ndarray) -> np.ndarray:
    """
    Compute backbone phi/psi dihedral angles.

    Args:
        coords: [N_res, 37, 3] float32

    Returns:
        [N_res, 4] float32 — sin_phi, cos_phi, sin_psi, cos_psi
        Undefined angles are set to 0.
    """
    n_res = coords.shape[0]
    phi_psi = np.zeros((n_res, 4), dtype=np.float32)

    N = coords[:, N_IDX, :]   # [N, 3]
    CA = coords[:, CA_IDX, :] # [N, 3]
    C = coords[:, C_IDX, :]   # [N, 3]

    def dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
        b1 = b - a
        b2 = c - b
        b3 = d - c
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-8))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        return float(np.arctan2(y, x))

    for i in range(n_res):
        # Phi: C(i-1) - N(i) - CA(i) - C(i)
        if i > 0 and not (
            np.any(np.isnan(C[i - 1])) or np.any(np.isnan(N[i]))
            or np.any(np.isnan(CA[i])) or np.any(np.isnan(C[i]))
        ):
            phi = dihedral(C[i - 1], N[i], CA[i], C[i])
            phi_psi[i, 0] = np.sin(phi)
            phi_psi[i, 1] = np.cos(phi)

        # Psi: N(i) - CA(i) - C(i) - N(i+1)
        if i < n_res - 1 and not (
            np.any(np.isnan(N[i])) or np.any(np.isnan(CA[i]))
            or np.any(np.isnan(C[i])) or np.any(np.isnan(N[i + 1]))
        ):
            psi = dihedral(N[i], CA[i], C[i], N[i + 1])
            phi_psi[i, 2] = np.sin(psi)
            phi_psi[i, 3] = np.cos(psi)

    return phi_psi


def run_dssp(pdb_path: Path) -> Optional[dict[int, str]]:
    """
    Run DSSP on a PDB file and return residue_index → SS_code mapping.

    Falls back to None if DSSP is unavailable.
    """
    try:
        from Bio.PDB import PDBParser, DSSP as BioDSSP
    except ImportError:
        logger.warning("BioPython DSSP not available")
        return None

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("prot", str(pdb_path))
        model = structure[0]
        dssp = BioDSSP(model, str(pdb_path), dssp="mkdssp")
        result: dict[int, str] = {}
        for key in dssp.keys():
            chain_id, (het, seq_id, icode) = key
            ss = dssp[key][2]
            result[seq_id] = ss if ss in DSSP_SS_MAP else " "
        return result
    except Exception as e:
        logger.debug(f"DSSP failed: {e}")
        return None


def dssp_to_sse_ids(ss_sequence: str) -> np.ndarray:
    """Convert a string of DSSP codes to integer SSE IDs [N_res] int8."""
    return np.array([DSSP_SS_MAP.get(c, DSSP_SS_MAP[" "]) for c in ss_sequence], dtype=np.int8)


def ss_to_simplified(ss_code: str) -> str:
    """Map full DSSP code to simplified H/E/C."""
    if ss_code in ("H", "G", "I"):
        return "H"
    if ss_code in ("E", "B"):
        return "E"
    return "C"
