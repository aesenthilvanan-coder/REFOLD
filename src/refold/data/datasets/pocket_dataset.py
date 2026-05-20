"""PyTorch dataset for pocket detector training."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


class PocketDataset:
    """
    Dataset of protein conformations with annotated binding pockets.

    Each item contains:
      - node_feats: [N_res, node_feat_dim] float32
      - edge_index: [2, n_edges] int64
      - edge_feats: [n_edges, edge_feat_dim] float32
      - pocket_labels: [N_res] float32 (1.0 = in pocket, 0.0 = not)
    """

    def __init__(
        self,
        structures: list,
        pocket_labels: list[dict],
        cache_dir: Optional[Path] = None,
    ):
        self.structures = structures
        self.pocket_labels = pocket_labels
        self.cache_dir = cache_dir
        self._cache: dict[int, dict] = {}

    def __len__(self) -> int:
        return len(self.structures)

    def __getitem__(self, idx: int) -> dict:
        if idx in self._cache:
            return self._cache[idx]

        structure = self.structures[idx]
        labels = self.pocket_labels[idx]

        try:
            item = self._build_item(structure, labels)
        except Exception as e:
            logger.warning(f"PocketDataset item {idx} failed: {e}")
            item = self._dummy_item()

        if self.cache_dir is not None:
            self._cache[idx] = item

        return item

    def _build_item(self, structure, labels: dict) -> dict:
        from refold.constants import GNN_NODE_DIM, GNN_EDGE_DIM, CONTACT_CA_CUTOFF
        from refold.constants import AA_TO_IDX

        n = structure.n_residues
        ca = structure.ca_coords

        # Node features
        node_feats = np.zeros((n, GNN_NODE_DIM), dtype=np.float32)
        for i in range(n):
            aa = structure.sequence[i] if i < len(structure.sequence) else "G"
            if aa in AA_TO_IDX:
                node_feats[i, AA_TO_IDX[aa]] = 1.0
            node_feats[i, 21] = float(structure.bfactors[i]) / 100.0
            node_feats[i, 22] = float(i) / max(n, 1)
            if structure.phi_psi is not None:
                node_feats[i, 23:27] = structure.phi_psi[i]
            if structure.rel_asa is not None:
                v = float(structure.rel_asa[i])
                node_feats[i, 27] = 0.0 if np.isnan(v) else v

        # Pocket labels per residue
        pocket_labels_arr = np.zeros(n, dtype=np.float32)
        for res_idx in labels.get("residue_indices", []):
            if 0 <= res_idx < n:
                pocket_labels_arr[res_idx] = 1.0

        # Edges from CA distance
        valid = ~np.any(np.isnan(ca), axis=-1)
        edges_i, edges_j, edge_feats_list = [], [], []

        for i in range(n):
            if not valid[i]:
                continue
            for j in range(i + 1, n):
                if not valid[j]:
                    continue
                d = float(np.linalg.norm(ca[i] - ca[j]))
                if d <= CONTACT_CA_CUTOFF:
                    ef = np.zeros(GNN_EDGE_DIM, dtype=np.float32)
                    ef[0] = d / 10.0
                    ef[1] = 1.0 / (d + 1e-8)
                    for k, center in enumerate(np.linspace(2.0, 12.0, 8)):
                        ef[2 + k] = float(np.exp(-((d - center) ** 2) / 4.0))
                    edges_i += [i, j]
                    edges_j += [j, i]
                    edge_feats_list += [ef, ef]

        if edges_i:
            edge_index = np.array([edges_i, edges_j], dtype=np.int64)
            edge_feats = np.array(edge_feats_list, dtype=np.float32)
        else:
            edge_index = np.zeros((2, 0), dtype=np.int64)
            edge_feats = np.zeros((0, GNN_EDGE_DIM), dtype=np.float32)

        return {
            "node_feats": node_feats,
            "edge_index": edge_index,
            "edge_feats": edge_feats,
            "pocket_labels": pocket_labels_arr,
        }

    def _dummy_item(self) -> dict:
        from refold.constants import GNN_NODE_DIM, GNN_EDGE_DIM
        return {
            "node_feats": np.zeros((10, GNN_NODE_DIM), dtype=np.float32),
            "edge_index": np.zeros((2, 0), dtype=np.int64),
            "edge_feats": np.zeros((0, GNN_EDGE_DIM), dtype=np.float32),
            "pocket_labels": np.zeros(10, dtype=np.float32),
        }
