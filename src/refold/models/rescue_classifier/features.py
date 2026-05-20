"""Feature extraction for the rescue classifier."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from refold.constants import (
    STANDARD_AAS, AA_TO_IDX, AA_PROPERTIES, CONTACT_CA_CUTOFF,
    THERMO_FEAT_DIM, EVO_FEAT_DIM, ESM2_EMBED_DIM, CA_IDX,
)
from refold.types import ProteinStructure, MutantStructure, Mutation

logger = logging.getLogger(__name__)


def build_thermodynamic_features(
    mutation: Mutation,
    ddg_predicted: float,
    ddg_esm1v: Optional[float] = None,
    ddg_gnn: Optional[float] = None,
) -> np.ndarray:
    """
    Build 16-dim thermodynamic feature vector.

    Features:
      [0]  ddg (clipped ±10)
      [1]  |ddg|
      [2]  ddg / 5.0 (normalized)
      [3]  ddg_esm1v or 0
      [4]  ddg_gnn or 0
      [5]  is_destabilizing (ddg > 0.5)
      [6]  is_severely_unstable (ddg > 5.0)
      [7]  wt_hydrophobicity
      [8]  mt_hydrophobicity
      [9]  delta_hydrophobicity
      [10] wt_charge
      [11] mt_charge
      [12] delta_charge
      [13] wt_volume
      [14] mt_volume
      [15] delta_volume_normalized
    """
    ddg = float(np.clip(ddg_predicted, -10, 10))
    ddg_e = float(ddg_esm1v) if ddg_esm1v is not None else 0.0
    ddg_g = float(ddg_gnn) if ddg_gnn is not None else 0.0

    wt = mutation.wildtype_aa
    mt = mutation.mutant_aa
    wt_props = AA_PROPERTIES.get(wt, {"hydrophobicity": 0.0, "charge": 0, "mw": 130.0, "volume": 130.0})
    mt_props = AA_PROPERTIES.get(mt, {"hydrophobicity": 0.0, "charge": 0, "mw": 130.0, "volume": 130.0})

    wt_h = wt_props["hydrophobicity"] / 5.0
    mt_h = mt_props["hydrophobicity"] / 5.0
    wt_v = wt_props["volume"] / 200.0
    mt_v = mt_props["volume"] / 200.0

    feat = np.array([
        ddg / 10.0,
        abs(ddg) / 10.0,
        ddg / 5.0,
        ddg_e / 10.0,
        ddg_g / 10.0,
        float(ddg > 0.5),
        float(ddg > 5.0),
        wt_h,
        mt_h,
        mt_h - wt_h,
        float(wt_props["charge"]),
        float(mt_props["charge"]),
        float(mt_props["charge"] - wt_props["charge"]),
        wt_v,
        mt_v,
        mt_v - wt_v,
    ], dtype=np.float32)

    assert len(feat) == THERMO_FEAT_DIM, f"Expected {THERMO_FEAT_DIM} thermo features, got {len(feat)}"
    return feat


def build_evolutionary_features(
    mutation: Mutation,
    structure: Optional[ProteinStructure] = None,
) -> np.ndarray:
    """
    Build 32-dim evolutionary/positional feature vector.

    Features derived from sequence position and amino acid properties.
    In a full implementation these would include MSA conservation scores.
    """
    n_aa = len(STANDARD_AAS)

    wt_onehot = np.zeros(n_aa, dtype=np.float32)
    mt_onehot = np.zeros(n_aa, dtype=np.float32)
    if mutation.wildtype_aa in AA_TO_IDX:
        wt_onehot[AA_TO_IDX[mutation.wildtype_aa]] = 1.0
    if mutation.mutant_aa in AA_TO_IDX:
        mt_onehot[AA_TO_IDX[mutation.mutant_aa]] = 1.0

    # Positional encoding (normalized)
    seq_len = structure.n_residues if structure is not None else 500
    pos_norm = float(mutation.position) / max(seq_len, 1)

    # Local secondary structure context (simplified)
    ss_feat = np.zeros(8, dtype=np.float32)
    if structure is not None and structure.sse_ids is not None:
        mut_idx = mutation.position - 1  # 0-based
        window = 4
        for k in range(-window, window + 1):
            idx = mut_idx + k
            if 0 <= idx < structure.n_residues:
                ss_code = int(structure.sse_ids[idx])
                if ss_code < 8:
                    ss_feat[ss_code] += 1
        ss_feat /= (ss_feat.sum() + 1e-8)

    # pLDDT at mutation site
    plddt_feat = 0.5
    if structure is not None:
        mut_idx = mutation.position - 1
        if 0 <= mut_idx < structure.n_residues:
            plddt_feat = float(structure.bfactors[mut_idx]) / 100.0

    # Relative solvent accessibility at mutation site
    rsa_feat = 0.5
    if structure is not None and structure.rel_asa is not None:
        mut_idx = mutation.position - 1
        if 0 <= mut_idx < structure.n_residues:
            v = float(structure.rel_asa[mut_idx])
            rsa_feat = 0.0 if np.isnan(v) else v

    feat = np.concatenate([
        wt_onehot[:6],   # partial wt onehot (first 6 properties)
        mt_onehot[:6],   # partial mt onehot
        ss_feat,          # 8 dim
        [pos_norm, plddt_feat, rsa_feat, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        np.zeros(3, dtype=np.float32),  # reserved for MSA features
    ], dtype=np.float32)

    # Pad or truncate to EVO_FEAT_DIM
    if len(feat) < EVO_FEAT_DIM:
        feat = np.pad(feat, (0, EVO_FEAT_DIM - len(feat)))
    else:
        feat = feat[:EVO_FEAT_DIM]

    assert len(feat) == EVO_FEAT_DIM, f"Expected {EVO_FEAT_DIM} evo features, got {len(feat)}"
    return feat


def build_graph_features(
    structure: ProteinStructure,
    mutation: Mutation,
    radius: float = CONTACT_CA_CUTOFF,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build local residue graph around mutation site.

    Returns:
        node_feats: [n_nodes, node_feat_dim] float32
        edge_feats: [n_edges, edge_feat_dim] float32
        edge_index:  [2, n_edges] int64
    """
    from refold.constants import GNN_NODE_DIM, GNN_EDGE_DIM

    mut_idx = mutation.position - 1
    n = structure.n_residues
    ca = structure.ca_coords  # [N, 3]

    # Select neighborhood
    mut_ca = ca[mut_idx] if not np.any(np.isnan(ca[mut_idx])) else np.zeros(3)
    dist_to_mut = np.sqrt(np.sum((ca - mut_ca) ** 2, axis=-1))
    in_radius = (dist_to_mut <= radius) & (~np.any(np.isnan(ca), axis=-1))
    node_indices = np.where(in_radius)[0]

    if len(node_indices) == 0:
        node_indices = np.array([mut_idx], dtype=np.int64)

    n_nodes = len(node_indices)
    node_feats = np.zeros((n_nodes, GNN_NODE_DIM), dtype=np.float32)

    for local_i, global_i in enumerate(node_indices):
        aa = structure.sequence[global_i] if global_i < len(structure.sequence) else "G"
        # One-hot AA (20 dims)
        if aa in AA_TO_IDX:
            node_feats[local_i, AA_TO_IDX[aa]] = 1.0
        # [20] is mutation site flag
        node_feats[local_i, 20] = float(global_i == mut_idx)
        # [21] pLDDT normalized
        node_feats[local_i, 21] = float(structure.bfactors[global_i]) / 100.0
        # [22] relative position
        node_feats[local_i, 22] = float(global_i) / max(n, 1)
        # [23-26] sin/cos dihedral if available
        if structure.phi_psi is not None and global_i < len(structure.phi_psi):
            node_feats[local_i, 23:27] = structure.phi_psi[global_i]
        # [27] rel_asa if available
        if structure.rel_asa is not None and global_i < len(structure.rel_asa):
            v = float(structure.rel_asa[global_i])
            node_feats[local_i, 27] = 0.0 if np.isnan(v) else v
        # [28] SSE type if available
        if structure.sse_ids is not None and global_i < len(structure.sse_ids):
            node_feats[local_i, 28] = float(structure.sse_ids[global_i]) / 7.0
        # [29] reserved
        node_feats[local_i, 29] = 0.0

    # Build edges (within-neighborhood contacts)
    local_ca = ca[node_indices]
    edges_i, edges_j = [], []
    for li in range(n_nodes):
        for lj in range(n_nodes):
            if li == lj:
                continue
            if np.any(np.isnan(local_ca[li])) or np.any(np.isnan(local_ca[lj])):
                continue
            d = float(np.linalg.norm(local_ca[li] - local_ca[lj]))
            if d <= radius:
                edges_i.append(li)
                edges_j.append(lj)

    if not edges_i:
        edge_index = np.zeros((2, 0), dtype=np.int64)
        edge_feats = np.zeros((0, GNN_EDGE_DIM), dtype=np.float32)
    else:
        edge_index = np.array([edges_i, edges_j], dtype=np.int64)
        n_edges = len(edges_i)
        edge_feats = np.zeros((n_edges, GNN_EDGE_DIM), dtype=np.float32)
        for e, (li, lj) in enumerate(zip(edges_i, edges_j)):
            gi, gj = node_indices[li], node_indices[lj]
            d = float(np.linalg.norm(ca[gi] - ca[gj]))
            # [0] distance / 10
            edge_feats[e, 0] = d / 10.0
            # [1] 1/d
            edge_feats[e, 1] = 1.0 / (d + 1e-8)
            # [2] Gaussian RBF
            for k, center in enumerate(np.linspace(2.0, 12.0, 8)):
                edge_feats[e, 2 + k] = float(np.exp(-((d - center) ** 2) / 4.0))
            # [10] sequence distance normalized
            edge_feats[e, 10] = float(abs(int(gi) - int(gj))) / max(n, 1)
            # [11] same SSE flag
            if structure.sse_ids is not None:
                edge_feats[e, 11] = float(structure.sse_ids[gi] == structure.sse_ids[gj])

    return node_feats, edge_feats, edge_index
