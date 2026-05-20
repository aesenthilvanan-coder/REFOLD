"""Rescue classifier model.

GNN (30-dim node, 20-dim edge) + ESM-2 (480-dim) +
thermodynamic (16-dim) + evolutionary (32-dim) → rescue probability.

Uses focal loss for 1:4 positive:negative class imbalance.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from refold.constants import (
    GNN_NODE_DIM, GNN_EDGE_DIM,
    ESM2_EMBED_DIM, THERMO_FEAT_DIM, EVO_FEAT_DIM,
    RESCUE_CLASSIFIER_HIDDEN_DIMS, RESCUE_CLASSIFIER_DROPOUT,
    FOCAL_LOSS_ALPHA, FOCAL_LOSS_GAMMA,
)


class GNNLayer(nn.Module):
    """Single message-passing layer with edge features."""

    def __init__(self, node_dim: int, edge_dim: int, hidden_dim: int):
        super().__init__()
        self.message_mlp = nn.Sequential(
            nn.Linear(node_dim * 2 + edge_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.update_mlp = nn.Sequential(
            nn.Linear(node_dim + hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
        )
    def forward(
        self,
        node_feats: torch.Tensor,
        edge_index: torch.Tensor,
        edge_feats: torch.Tensor,
    ) -> torch.Tensor:
        src, dst = edge_index[0], edge_index[1]
        msg_input = torch.cat([node_feats[src], node_feats[dst], edge_feats], dim=-1)
        messages = self.message_mlp(msg_input)

        agg = torch.zeros(
            node_feats.shape[0], messages.shape[-1],
            device=node_feats.device, dtype=node_feats.dtype
        )
        agg.scatter_add_(0, dst.unsqueeze(-1).expand_as(messages), messages)

        updated = self.update_mlp(torch.cat([node_feats, agg], dim=-1))
        return updated


class StructureGNN(nn.Module):
    """Multi-layer GNN for protein structure encoding."""

    def __init__(
        self,
        node_in_dim: int = GNN_NODE_DIM,
        edge_in_dim: int = GNN_EDGE_DIM,
        hidden_dim: int = 256,
        n_layers: int = 4,
        output_dim: int = 256,
    ):
        super().__init__()
        self.layers = nn.ModuleList()
        in_dim = node_in_dim
        for _ in range(n_layers):
            self.layers.append(GNNLayer(in_dim, edge_in_dim, hidden_dim))
            in_dim = hidden_dim

        self.pool = nn.Linear(hidden_dim, output_dim)

    def forward(
        self,
        node_feats: torch.Tensor,
        edge_index: torch.Tensor,
        edge_feats: torch.Tensor,
        node_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        # node_feats: [B, N_max, node_dim] or [N_total, node_dim]
        if node_feats.dim() == 3:
            # Batched: flatten to [N_total, dim]
            B, N, D = node_feats.shape
            x = node_feats.view(B * N, D)
            for layer in self.layers:
                x = layer(x, edge_index, edge_feats)
            x = x.view(B, N, -1)
            if node_mask is not None:
                x = x * node_mask.unsqueeze(-1).float()
                denom = node_mask.float().sum(dim=1, keepdim=True).clamp(min=1)
                graph_repr = x.sum(dim=1) / denom
            else:
                graph_repr = x.mean(dim=1)
        else:
            x = node_feats
            for layer in self.layers:
                x = layer(x, edge_index, edge_feats)
            graph_repr = x.mean(dim=0, keepdim=True)

        return self.pool(graph_repr)


class RescueClassifier(nn.Module):
    """Full rescue amenability classifier.

    Fuses structural GNN + ESM-2 embedding + thermodynamic + evolutionary features.
    """

    def __init__(
        self,
        gnn_node_dim: int = GNN_NODE_DIM,
        gnn_edge_dim: int = GNN_EDGE_DIM,
        gnn_hidden_dim: int = 256,
        gnn_n_layers: int = 4,
        gnn_output_dim: int = 256,
        esm2_dim: int = ESM2_EMBED_DIM,
        thermo_dim: int = THERMO_FEAT_DIM,
        evo_dim: int = EVO_FEAT_DIM,
        hidden_dims: list[int] = RESCUE_CLASSIFIER_HIDDEN_DIMS,
        dropout: float = RESCUE_CLASSIFIER_DROPOUT,
    ):
        super().__init__()

        self.structure_gnn = StructureGNN(
            node_in_dim=gnn_node_dim,
            edge_in_dim=gnn_edge_dim,
            hidden_dim=gnn_hidden_dim,
            n_layers=gnn_n_layers,
            output_dim=gnn_output_dim,
        )

        self.esm2_proj = nn.Sequential(
            nn.Linear(esm2_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
        )

        self.thermo_proj = nn.Sequential(
            nn.Linear(thermo_dim, 64),
            nn.GELU(),
        )

        self.evo_proj = nn.Sequential(
            nn.Linear(evo_dim, 64),
            nn.GELU(),
        )

        fusion_dim = gnn_output_dim + 256 + 64 + 64
        layers = []
        in_dim = fusion_dim
        for out_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_dim, out_dim),
                nn.LayerNorm(out_dim),
                nn.GELU(),
                nn.Dropout(dropout),
            ])
            in_dim = out_dim

        layers.append(nn.Linear(in_dim, 1))
        self.classifier = nn.Sequential(*layers)

    def forward(
        self,
        node_feats: torch.Tensor,
        edge_index: torch.Tensor,
        edge_feats: torch.Tensor,
        esm2_embedding: torch.Tensor,
        thermo_features: torch.Tensor,
        evo_features: torch.Tensor,
        node_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward pass.

        Returns logits [B, 1]. Apply sigmoid for probability.
        """
        gnn_repr = self.structure_gnn(node_feats, edge_index, edge_feats, node_mask)
        esm2_repr = self.esm2_proj(esm2_embedding)
        thermo_repr = self.thermo_proj(thermo_features)
        evo_repr = self.evo_proj(evo_features)

        if gnn_repr.shape[0] != esm2_repr.shape[0]:
            gnn_repr = gnn_repr.expand(esm2_repr.shape[0], -1)

        fusion = torch.cat([gnn_repr, esm2_repr, thermo_repr, evo_repr], dim=-1)
        logits = self.classifier(fusion)
        return logits

    def predict_proba(
        self,
        node_feats: torch.Tensor,
        edge_index: torch.Tensor,
        edge_feats: torch.Tensor,
        esm2_embedding: torch.Tensor,
        thermo_features: torch.Tensor,
        evo_features: torch.Tensor,
    ) -> torch.Tensor:
        logits = self(
            node_feats, edge_index, edge_feats,
            esm2_embedding, thermo_features, evo_features,
        )
        return torch.sigmoid(logits)


class FocalLoss(nn.Module):
    """Focal loss for handling class imbalance."""

    def __init__(
        self,
        alpha: float = FOCAL_LOSS_ALPHA,
        gamma: float = FOCAL_LOSS_GAMMA,
        reduction: str = "mean",
    ):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        bce_loss = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
        probs = torch.sigmoid(logits)
        p_t = probs * targets + (1 - probs) * (1 - targets)
        alpha_t = self.alpha * targets + (1 - self.alpha) * (1 - targets)
        focal_weight = alpha_t * (1 - p_t) ** self.gamma
        loss = focal_weight * bce_loss

        if self.reduction == "mean":
            return loss.mean()
        elif self.reduction == "sum":
            return loss.sum()
        return loss
