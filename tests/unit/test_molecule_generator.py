"""Unit tests for molecule generator diffusion model."""

import torch
import pytest

from refold.models.molecule_generator.diffusion import (
    NoiseSchedule, PocketEncoder, REFOLDDiffusionModel,
)
from refold.constants import (
    N_MOL_ATOM_TYPES, MOL_GEN_POCKET_EMBED_DIM,
    POCKET_NODE_FEAT_DIM, POCKET_SCALAR_DIM,
)


def test_noise_schedule_q_sample():
    schedule = NoiseSchedule(T=100)
    x0 = torch.randn(2, 10, 3)
    t = torch.tensor([0, 50])
    xt = schedule.q_sample(x0, t)
    assert xt.shape == x0.shape
    assert not torch.allclose(xt, x0)


def test_noise_schedule_p_sample():
    schedule = NoiseSchedule(T=100)
    x_t = torch.randn(2, 10, 3)
    model_out = torch.randn(2, 10, 3)
    x_prev = schedule.p_sample(model_out, x_t, t=50)
    assert x_prev.shape == x_t.shape


def test_pocket_encoder_forward():
    encoder = PocketEncoder(output_dim=MOL_GEN_POCKET_EMBED_DIM)
    n_spheres = 20
    sphere_feats = torch.randn(n_spheres, POCKET_NODE_FEAT_DIM)
    scalars = torch.randn(1, POCKET_SCALAR_DIM)
    emb = encoder(sphere_feats, scalars)
    assert emb.shape == (MOL_GEN_POCKET_EMBED_DIM,)


def test_diffusion_model_forward(dummy_pocket):
    model = REFOLDDiffusionModel(T=10)
    model.eval()
    B, N_atoms = 2, 8
    coords = torch.randn(B, N_atoms, 3)
    atom_types = torch.zeros(B, N_atoms, N_MOL_ATOM_TYPES)
    atom_types[:, :, 0] = 1.0
    from refold.models.molecule_generator.diffusion import _build_pocket_features
    sphere_feats, scalars = _build_pocket_features(dummy_pocket)
    sphere_feats = sphere_feats.unsqueeze(0).expand(B, -1, -1)
    scalars_b = scalars.expand(B, -1)
    t = torch.randint(0, 10, (B,))
    with torch.no_grad():
        loss = model(coords, atom_types, sphere_feats, scalars_b, t)
    assert loss.item() >= 0
