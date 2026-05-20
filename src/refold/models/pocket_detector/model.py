"""Pocket detector neural network model."""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from refold.constants import (
    GNN_NODE_DIM, GNN_EDGE_DIM,
    POCKET_GNN_HIDDEN_DIM, POCKET_GNN_N_LAYERS, POCKET_GNN_N_HEADS,
)


class PocketDetectorGNN(nn.Module):
    """GNN for residue-level pocket membership prediction."""

    def __init__(
        self,
        node_in_dim: int = GNN_NODE_DIM,
        edge_in_dim: int = GNN_EDGE_DIM,
        hidden_dim: int = POCKET_GNN_HIDDEN_DIM,
        n_layers: int = POCKET_GNN_N_LAYERS,
        n_heads: int = POCKET_GNN_N_HEADS,
    ):
        super().__init__()
        self.node_proj = nn.Linear(node_in_dim, hidden_dim)
        self.edge_proj = nn.Linear(edge_in_dim, hidden_dim)

        self.attention_layers = nn.ModuleList([
            nn.MultiheadAttention(
                embed_dim=hidden_dim,
                num_heads=n_heads,
                batch_first=True,
                dropout=0.1,
            )
            for _ in range(n_layers)
        ])
        self.ff_layers = nn.ModuleList([
            nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim * 4),
                nn.GELU(),
                nn.Linear(hidden_dim * 4, hidden_dim),
                nn.LayerNorm(hidden_dim),
            )
            for _ in range(n_layers)
        ])
        self.layer_norms = nn.ModuleList([
            nn.LayerNorm(hidden_dim) for _ in range(n_layers)
        ])

        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, 64),
            nn.GELU(),
            nn.Linear(64, 1),
        )

    def forward(
        self,
        node_feats: torch.Tensor,
        edge_index: torch.Tensor,
        edge_feats: torch.Tensor,
        node_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward pass for pocket residue classification.

        Args:
            node_feats: [B, N, node_in_dim] or [N, node_in_dim]
            Returns: [B, N, 1] or [N, 1] pocket membership logits.
        """
        if node_feats.dim() == 2:
            node_feats = node_feats.unsqueeze(0)
            squeeze = True
        else:
            squeeze = False

        x = self.node_proj(node_feats)

        key_padding_mask = None
        if node_mask is not None:
            key_padding_mask = ~node_mask

        for attn, ff, ln in zip(self.attention_layers, self.ff_layers, self.layer_norms):
            attn_out, _ = attn(x, x, x, key_padding_mask=key_padding_mask)
            x = ln(x + attn_out)
            x = x + ff(x)

        logits = self.classifier(x)
        if squeeze:
            logits = logits.squeeze(0)
        return logits
