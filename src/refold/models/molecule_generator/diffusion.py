"""Pocket-conditioned DDPM diffusion model for de novo molecule generation."""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from refold.constants import (
    N_DIFFUSION_STEPS_DEFAULT,
    DIFFUSION_BETA_START,
    DIFFUSION_BETA_END,
    N_MOL_ATOM_TYPES,
    MOL_ATOM_TYPES,
    MOL_GEN_HIDDEN_DIM,
    MOL_GEN_N_LAYERS,
    MOL_GEN_POCKET_EMBED_DIM,
    POCKET_NODE_FEAT_DIM,
    POCKET_SCALAR_DIM,
    N_ATOMS_MIN,
    N_ATOMS_MAX,
)

if TYPE_CHECKING:
    from refold.types import Pocket

logger = logging.getLogger(__name__)


class NoiseSchedule(nn.Module):
    """Linear beta noise schedule for DDPM."""

    def __init__(self, T: int = N_DIFFUSION_STEPS_DEFAULT):
        super().__init__()
        self.T = T
        betas = torch.linspace(DIFFUSION_BETA_START, DIFFUSION_BETA_END, T)
        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        alphas_cumprod_prev = F.pad(alphas_cumprod[:-1], (1, 0), value=1.0)

        self.register_buffer("betas", betas)
        self.register_buffer("alphas", alphas)
        self.register_buffer("alphas_cumprod", alphas_cumprod)
        self.register_buffer("alphas_cumprod_prev", alphas_cumprod_prev)
        self.register_buffer("sqrt_alphas_cumprod", torch.sqrt(alphas_cumprod))
        self.register_buffer("sqrt_one_minus_alphas_cumprod", torch.sqrt(1.0 - alphas_cumprod))
        self.register_buffer(
            "posterior_variance",
            betas * (1.0 - alphas_cumprod_prev) / (1.0 - alphas_cumprod),
        )

    def q_sample(
        self,
        x_0: torch.Tensor,
        t: torch.Tensor,
        noise: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward diffusion: sample x_t from x_0 at timestep t."""
        if noise is None:
            noise = torch.randn_like(x_0)
        sqrt_alpha = self.sqrt_alphas_cumprod[t].view(-1, *([1] * (x_0.dim() - 1)))
        sqrt_one_minus = self.sqrt_one_minus_alphas_cumprod[t].view(-1, *([1] * (x_0.dim() - 1)))
        return sqrt_alpha * x_0 + sqrt_one_minus * noise

    @torch.no_grad()
    def p_sample(
        self,
        model_out: torch.Tensor,
        x_t: torch.Tensor,
        t: int,
    ) -> torch.Tensor:
        """Single reverse diffusion step: sample x_{t-1} from x_t."""
        beta_t = self.betas[t]
        sqrt_recip_alpha = 1.0 / torch.sqrt(self.alphas[t])
        sqrt_one_minus = self.sqrt_one_minus_alphas_cumprod[t]

        mean = sqrt_recip_alpha * (x_t - beta_t / sqrt_one_minus * model_out)

        if t == 0:
            return mean
        else:
            noise = torch.randn_like(x_t)
            var = self.posterior_variance[t]
            return mean + torch.sqrt(var) * noise


class PocketEncoder(nn.Module):
    """Encode pocket geometry into a fixed-dim embedding."""

    def __init__(
        self,
        node_feat_dim: int = POCKET_NODE_FEAT_DIM,
        scalar_dim: int = POCKET_SCALAR_DIM,
        output_dim: int = MOL_GEN_POCKET_EMBED_DIM,
        hidden_dim: int = 128,
    ):
        super().__init__()

        self.sphere_mlp = nn.Sequential(
            nn.Linear(node_feat_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        self.scalar_mlp = nn.Sequential(
            nn.Linear(scalar_dim, 32),
            nn.GELU(),
        )

        self.output_proj = nn.Sequential(
            nn.Linear(hidden_dim + 32, output_dim),
            nn.LayerNorm(output_dim),
            nn.GELU(),
        )

    def forward(
        self,
        alpha_sphere_feats: torch.Tensor,
        pocket_scalars: torch.Tensor,
    ) -> torch.Tensor:
        """Encode pocket.

        Args:
            alpha_sphere_feats: [n_spheres, node_feat_dim] or [B, n_spheres, node_feat_dim]
            pocket_scalars: [1, scalar_dim] or [B, scalar_dim]

        Returns: [output_dim] or [B, output_dim]
        """
        if alpha_sphere_feats.dim() == 2:
            sphere_enc = self.sphere_mlp(alpha_sphere_feats).mean(dim=0, keepdim=True)
            scalar_enc = self.scalar_mlp(pocket_scalars)
            combined = torch.cat([sphere_enc, scalar_enc], dim=-1)
            return self.output_proj(combined).squeeze(0)
        else:
            sphere_enc = self.sphere_mlp(alpha_sphere_feats).mean(dim=1)
            scalar_enc = self.scalar_mlp(pocket_scalars)
            combined = torch.cat([sphere_enc, scalar_enc], dim=-1)
            return self.output_proj(combined)


class DenoisingNetwork(nn.Module):
    """Denoising network for atom coordinate and type prediction."""

    def __init__(
        self,
        n_atom_types: int = N_MOL_ATOM_TYPES,
        hidden_dim: int = MOL_GEN_HIDDEN_DIM,
        n_layers: int = MOL_GEN_N_LAYERS,
        pocket_embed_dim: int = MOL_GEN_POCKET_EMBED_DIM,
        n_atoms_max: int = N_ATOMS_MAX,
    ):
        super().__init__()
        self.n_atom_types = n_atom_types
        self.n_atoms_max = n_atoms_max

        # Atom feature encoder: 3 (coords) + n_atom_types (one-hot) + 1 (time embedding)
        atom_in_dim = 3 + n_atom_types

        # Time embedding: 1 → hidden_dim
        self.time_emb = nn.Sequential(
            nn.Linear(1, hidden_dim // 4),
            nn.GELU(),
            nn.Linear(hidden_dim // 4, hidden_dim),
        )

        self.atom_proj = nn.Linear(atom_in_dim, hidden_dim)
        self.pocket_proj = nn.Linear(pocket_embed_dim, hidden_dim)

        # Transformer layers for atom interactions
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_dim,
            nhead=8,
            dim_feedforward=hidden_dim * 4,
            dropout=0.1,
            batch_first=True,
            norm_first=True,
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        # Output heads
        self.coords_head = nn.Linear(hidden_dim, 3)
        self.type_head = nn.Linear(hidden_dim, n_atom_types)

    def forward(
        self,
        coords: torch.Tensor,
        atom_types: torch.Tensor,
        t: torch.Tensor,
        pocket_embedding: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Predict noise for coordinates and atom types.

        Args:
            coords: [B, N_atoms, 3]
            atom_types: [B, N_atoms, n_atom_types] one-hot
            t: [B] timestep
            pocket_embedding: [B, pocket_embed_dim]

        Returns:
            (noise_coords [B, N_atoms, 3], noise_types [B, N_atoms, n_atom_types])
        """
        B, N, _ = coords.shape

        t_norm = (t.float() / max(self.n_atoms_max, 1)).unsqueeze(-1).unsqueeze(-1)
        t_emb = self.time_emb(t_norm).squeeze(-2)

        atom_in = torch.cat([coords, atom_types], dim=-1)
        x = self.atom_proj(atom_in)

        p_emb = self.pocket_proj(pocket_embedding).unsqueeze(1)
        x = x + p_emb
        x = x + t_emb.unsqueeze(1)

        x = self.transformer(x)

        noise_coords = self.coords_head(x)
        noise_types = self.type_head(x)
        return noise_coords, noise_types


def _build_pocket_features(pocket: "Pocket") -> tuple[torch.Tensor, torch.Tensor]:
    """Build pocket features for the encoder."""
    from refold.types import PocketType

    alpha_spheres = pocket.alpha_sphere_coords
    n_spheres = len(alpha_spheres)

    # Node features: [n_spheres, POCKET_NODE_FEAT_DIM=25]
    center = pocket.center
    node_feats = np.zeros((n_spheres, POCKET_NODE_FEAT_DIM), dtype=np.float32)

    # Normalized coordinates relative to pocket center
    rel_coords = alpha_spheres - center
    node_feats[:, :3] = rel_coords / 10.0

    # Distance to center
    node_feats[:, 3] = np.linalg.norm(rel_coords, axis=-1) / 10.0

    # Global scalars appended per-sphere
    node_feats[:, 4] = pocket.volume / 1000.0
    node_feats[:, 5] = pocket.druggability_score
    node_feats[:, 6] = pocket.hydrophobicity
    node_feats[:, 7] = pocket.detection_frequency
    node_feats[:, 8] = float(pocket.is_transient)

    sphere_tensor = torch.tensor(node_feats)
    scalar_tensor = _build_pocket_scalars(pocket)
    return sphere_tensor, scalar_tensor


def _build_pocket_scalars(pocket: "Pocket") -> torch.Tensor:
    """Build [1, 8] global scalar features for the pocket."""
    center = pocket.center
    scalars = np.array([
        pocket.volume / 1000.0,
        pocket.druggability_score,
        pocket.hydrophobicity / 5.0,
        pocket.detection_frequency,
        center[0] / 50.0,
        center[1] / 50.0,
        center[2] / 50.0,
        float(pocket.is_transient),
    ], dtype=np.float32)
    return torch.tensor(scalars).unsqueeze(0)


class REFOLDDiffusionModel(nn.Module):
    """Full pocket-conditioned diffusion model for de novo molecule generation."""

    def __init__(
        self,
        T: int = N_DIFFUSION_STEPS_DEFAULT,
        n_atom_types: int = N_MOL_ATOM_TYPES,
        hidden_dim: int = MOL_GEN_HIDDEN_DIM,
        n_layers: int = MOL_GEN_N_LAYERS,
        pocket_embed_dim: int = MOL_GEN_POCKET_EMBED_DIM,
    ):
        super().__init__()
        self.T = T
        self.n_atom_types = n_atom_types

        self.noise_schedule = NoiseSchedule(T=T)
        self.pocket_encoder = PocketEncoder(output_dim=pocket_embed_dim)
        self.denoising_net = DenoisingNetwork(
            n_atom_types=n_atom_types,
            hidden_dim=hidden_dim,
            n_layers=n_layers,
            pocket_embed_dim=pocket_embed_dim,
        )

    def forward(
        self,
        coords: torch.Tensor,
        atom_types: torch.Tensor,
        pocket_sphere_feats: torch.Tensor,
        pocket_scalars: torch.Tensor,
        t: torch.Tensor,
    ) -> torch.Tensor:
        """Training forward pass: predict noise given noisy inputs at timestep t."""
        pocket_emb = self.pocket_encoder(pocket_sphere_feats, pocket_scalars)

        noisy_coords = self.noise_schedule.q_sample(coords, t)
        noise_coords_pred, noise_types_pred = self.denoising_net(
            noisy_coords, atom_types, t, pocket_emb
        )

        noise_coords_true = torch.randn_like(coords)
        loss_coords = F.mse_loss(noise_coords_pred, noise_coords_true)
        loss_types = F.cross_entropy(
            noise_types_pred.view(-1, self.n_atom_types),
            atom_types.argmax(dim=-1).view(-1),
        )
        return loss_coords + 0.1 * loss_types

    @torch.no_grad()
    def sample(
        self,
        pocket: "Pocket",
        n_molecules: int = 100,
        guidance_scale: float = 1.5,
        n_atoms_range: tuple[int, int] = (N_ATOMS_MIN, N_ATOMS_MAX),
        device: torch.device | None = None,
    ) -> list[tuple[np.ndarray, list[str]]]:
        """Sample molecules conditioned on a pocket.

        Returns list of (coords [N_atoms, 3], atom_type_list) tuples.
        """
        if device is None:
            device = next(self.parameters()).device

        sphere_feats, scalars = _build_pocket_features(pocket)
        sphere_feats = sphere_feats.to(device).unsqueeze(0)
        scalars = scalars.to(device)
        pocket_emb = self.pocket_encoder(sphere_feats.squeeze(0), scalars)

        results = []
        rng = np.random.default_rng()

        for _ in range(n_molecules):
            n_atoms = int(rng.integers(n_atoms_range[0], n_atoms_range[1] + 1))

            coords = torch.randn(1, n_atoms, 3, device=device)
            atom_types = torch.zeros(1, n_atoms, self.n_atom_types, device=device)
            atom_types[:, :, 0] = 1.0

            for t_step in reversed(range(self.T)):
                t_tensor = torch.tensor([t_step], device=device)
                noise_c, noise_t = self.denoising_net(
                    coords, atom_types, t_tensor, pocket_emb.unsqueeze(0)
                )
                coords = self.noise_schedule.p_sample(noise_c, coords, t_step)
                atom_types = F.softmax(
                    F.softmax(atom_types, dim=-1) - 0.1 * noise_t, dim=-1
                )

            atom_type_indices = atom_types.squeeze(0).argmax(dim=-1).cpu().numpy()
            atom_type_names = [MOL_ATOM_TYPES[i] for i in atom_type_indices]
            coords_np = coords.squeeze(0).cpu().numpy()

            results.append((coords_np, atom_type_names))

        return results
