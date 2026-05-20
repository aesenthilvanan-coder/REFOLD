"""PyTorch Dataset for rescue classifier training."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class RescueDataset:
    """Dataset for rescue classifier — loads features per mutation, caches as .pt files."""

    def __init__(
        self,
        df: pd.DataFrame,
        structure_dir: Path,
        cache_dir: Path,
        esm2_embed_dim: int = 480,
    ):
        self.df = df.reset_index(drop=True)
        self.structure_dir = structure_dir
        self.cache_dir = cache_dir
        self.esm2_embed_dim = esm2_embed_dim
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, idx: int) -> dict[str, Any]:
        import torch

        row = self.df.iloc[idx]
        cache_key = (
            f"{row['uniprot_id']}_{row['position']}_"
            f"{row['wildtype_aa']}_{row['mutant_aa']}.pt"
        )
        cache_path = self.cache_dir / cache_key

        if cache_path.exists():
            try:
                return torch.load(cache_path, map_location="cpu")
            except Exception:
                pass

        item = self._compute_item(row)
        try:
            torch.save(item, cache_path)
        except Exception as e:
            logger.debug(f"Failed to cache item {cache_key}: {e}")
        return item

    def _compute_item(self, row: pd.Series) -> dict[str, Any]:
        import torch

        try:
            from refold.structure.parser import parse_pdb_to_structure
            from refold.data.alphafold_fetcher import get_structure_path
            from refold.constants import (
                GNN_NODE_DIM, GNN_EDGE_DIM, ESM2_EMBED_DIM,
                THERMO_FEAT_DIM, EVO_FEAT_DIM, CONTACT_CA_CUTOFF,
            )

            pdb_path = get_structure_path(row["uniprot_id"])
            if pdb_path is None:
                return self._dummy_item(row)

            structure = parse_pdb_to_structure(pdb_path, row["uniprot_id"])
            pos_0 = int(row["position"]) - 1

            if pos_0 < 0 or pos_0 >= structure.n_residues:
                return self._dummy_item(row)

            # Build structural graph around mutation site
            ca = structure.ca_coords  # [N, 3]
            mut_ca = ca[pos_0]
            dists = np.linalg.norm(ca - mut_ca, axis=-1)
            neighbor_mask = (dists < CONTACT_CA_CUTOFF) & structure.residue_mask

            neighbor_indices = np.where(neighbor_mask)[0]
            n_nodes = len(neighbor_indices)

            if n_nodes == 0:
                return self._dummy_item(row)

            # Node features: [n_nodes, GNN_NODE_DIM]
            node_feats = np.zeros((n_nodes, GNN_NODE_DIM), dtype=np.float32)
            for i, res_idx in enumerate(neighbor_indices):
                res_type = structure.residue_types[res_idx]
                res_type_idx = min(res_type, 19)
                node_feats[i, res_type_idx] = 1.0
                node_feats[i, 20] = structure.bfactors[res_idx] / 100.0
                if structure.rel_asa is not None:
                    node_feats[i, 21] = structure.rel_asa[res_idx]
                node_feats[i, 22] = float(res_idx == pos_0)
                if structure.sse_ids is not None:
                    sse = min(int(structure.sse_ids[res_idx]), 7)
                    node_feats[i, 23 + sse] = 1.0

            # Edge index and features
            src_list, dst_list, edge_feats_list = [], [], []
            for i, ri in enumerate(neighbor_indices):
                for j, rj in enumerate(neighbor_indices):
                    if i == j:
                        continue
                    d = float(np.linalg.norm(ca[ri] - ca[rj]))
                    if d < CONTACT_CA_CUTOFF:
                        src_list.append(i)
                        dst_list.append(j)
                        ef = np.zeros(GNN_EDGE_DIM, dtype=np.float32)
                        ef[0] = d / CONTACT_CA_CUTOFF
                        ef[1] = float(abs(ri - rj)) / structure.n_residues
                        edge_feats_list.append(ef)

            if not src_list:
                src_list, dst_list = [0], [0]
                edge_feats_list = [np.zeros(GNN_EDGE_DIM, dtype=np.float32)]

            label = float(row.get("rescue_label", 0.0))
            ddg = float(row.get("ddg_kcal_mol", 0.0))

            thermo_feat = np.zeros(THERMO_FEAT_DIM, dtype=np.float32)
            thermo_feat[0] = np.clip(ddg / 10.0, -2.0, 2.0)

            evo_feat = np.zeros(EVO_FEAT_DIM, dtype=np.float32)
            esm2_emb = np.zeros(self.esm2_embed_dim, dtype=np.float32)

            return {
                "node_feats": torch.tensor(node_feats),
                "edge_index": torch.tensor([src_list, dst_list], dtype=torch.long),
                "edge_feats": torch.tensor(np.array(edge_feats_list)),
                "esm2_embedding": torch.tensor(esm2_emb),
                "thermo_features": torch.tensor(thermo_feat),
                "evo_features": torch.tensor(evo_feat),
                "label": torch.tensor(label, dtype=torch.float32),
                "ddg": torch.tensor(ddg, dtype=torch.float32),
                "uniprot_id": str(row["uniprot_id"]),
                "position": int(row["position"]),
            }
        except Exception as e:
            logger.debug(f"Failed to compute item: {e}")
            return self._dummy_item(row)

    def _dummy_item(self, row: pd.Series) -> dict[str, Any]:
        import torch
        from refold.constants import GNN_NODE_DIM, GNN_EDGE_DIM, THERMO_FEAT_DIM, EVO_FEAT_DIM

        return {
            "node_feats": torch.zeros((1, GNN_NODE_DIM)),
            "edge_index": torch.zeros((2, 1), dtype=torch.long),
            "edge_feats": torch.zeros((1, GNN_EDGE_DIM)),
            "esm2_embedding": torch.zeros(self.esm2_embed_dim),
            "thermo_features": torch.zeros(THERMO_FEAT_DIM),
            "evo_features": torch.zeros(EVO_FEAT_DIM),
            "label": torch.tensor(float(row.get("rescue_label", 0.0))),
            "ddg": torch.tensor(float(row.get("ddg_kcal_mol", 0.0))),
            "uniprot_id": str(row.get("uniprot_id", "")),
            "position": int(row.get("position", 0)),
        }


def collate_rescue_batch(batch: list[dict[str, Any]]) -> dict[str, Any]:
    """Collate rescue dataset items — pads node/edge tensors to same size."""
    import torch

    max_nodes = max(item["node_feats"].shape[0] for item in batch)
    node_dim = batch[0]["node_feats"].shape[1]
    edge_dim = batch[0]["edge_feats"].shape[1]
    esm2_dim = batch[0]["esm2_embedding"].shape[0]
    thermo_dim = batch[0]["thermo_features"].shape[0]
    evo_dim = batch[0]["evo_features"].shape[0]

    node_feats_padded = torch.zeros(len(batch), max_nodes, node_dim)
    node_masks = torch.zeros(len(batch), max_nodes, dtype=torch.bool)
    esm2_embs = torch.zeros(len(batch), esm2_dim)
    thermo_feats = torch.zeros(len(batch), thermo_dim)
    evo_feats = torch.zeros(len(batch), evo_dim)
    labels = torch.zeros(len(batch))
    ddgs = torch.zeros(len(batch))

    edge_indices_list = []
    edge_feats_list = []
    batch_offsets = []

    offset = 0
    for i, item in enumerate(batch):
        n = item["node_feats"].shape[0]
        node_feats_padded[i, :n] = item["node_feats"]
        node_masks[i, :n] = True
        esm2_embs[i] = item["esm2_embedding"]
        thermo_feats[i] = item["thermo_features"]
        evo_feats[i] = item["evo_features"]
        labels[i] = item["label"]
        ddgs[i] = item["ddg"]

        ei = item["edge_index"] + offset
        edge_indices_list.append(ei)
        edge_feats_list.append(item["edge_feats"])
        batch_offsets.append(offset)
        offset += n

    return {
        "node_feats": node_feats_padded,
        "node_mask": node_masks,
        "edge_index": torch.cat(edge_indices_list, dim=1),
        "edge_feats": torch.cat(edge_feats_list, dim=0),
        "esm2_embedding": esm2_embs,
        "thermo_features": thermo_feats,
        "evo_features": evo_feats,
        "labels": labels,
        "ddgs": ddgs,
    }
