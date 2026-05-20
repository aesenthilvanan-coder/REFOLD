"""Diffusion noise schedule for molecule generation."""

from __future__ import annotations

import numpy as np

from refold.constants import (
    N_DIFFUSION_STEPS_DEFAULT, N_DIFFUSION_STEPS_M1,
    DIFFUSION_BETA_START, DIFFUSION_BETA_END,
)


def linear_beta_schedule(
    T: int = N_DIFFUSION_STEPS_DEFAULT,
    beta_start: float = DIFFUSION_BETA_START,
    beta_end: float = DIFFUSION_BETA_END,
) -> np.ndarray:
    """Return [T] linear beta schedule."""
    return np.linspace(beta_start, beta_end, T, dtype=np.float64)


def cosine_beta_schedule(T: int = N_DIFFUSION_STEPS_DEFAULT, s: float = 0.008) -> np.ndarray:
    """Return [T] cosine beta schedule (Nichol & Dhariwal 2021)."""
    steps = np.arange(T + 1)
    f = np.cos((steps / T + s) / (1 + s) * np.pi / 2) ** 2
    alphas_cumprod = f / f[0]
    betas = 1 - alphas_cumprod[1:] / alphas_cumprod[:-1]
    return np.clip(betas, 1e-8, 0.999).astype(np.float64)


def precompute_schedule_tensors(betas: np.ndarray) -> dict:
    """
    Precompute all derived noise schedule tensors.

    Returns dict with numpy arrays:
        betas, alphas, alphas_cumprod, alphas_cumprod_prev,
        sqrt_alphas_cumprod, sqrt_one_minus_alphas_cumprod,
        posterior_variance
    """
    T = len(betas)
    alphas = 1.0 - betas
    alphas_cumprod = np.cumprod(alphas)
    alphas_cumprod_prev = np.concatenate([[1.0], alphas_cumprod[:-1]])

    sqrt_alphas_cumprod = np.sqrt(alphas_cumprod)
    sqrt_one_minus_alphas_cumprod = np.sqrt(1.0 - alphas_cumprod)
    sqrt_recip_alphas = 1.0 / np.sqrt(alphas)

    # q(x_{t-1} | x_t, x_0) variance
    posterior_variance = betas * (1.0 - alphas_cumprod_prev) / (1.0 - alphas_cumprod)

    return {
        "T": T,
        "betas": betas.astype(np.float32),
        "alphas": alphas.astype(np.float32),
        "alphas_cumprod": alphas_cumprod.astype(np.float32),
        "alphas_cumprod_prev": alphas_cumprod_prev.astype(np.float32),
        "sqrt_alphas_cumprod": sqrt_alphas_cumprod.astype(np.float32),
        "sqrt_one_minus_alphas_cumprod": sqrt_one_minus_alphas_cumprod.astype(np.float32),
        "sqrt_recip_alphas": sqrt_recip_alphas.astype(np.float32),
        "posterior_variance": posterior_variance.astype(np.float32),
    }


def get_default_schedule(T: int = N_DIFFUSION_STEPS_DEFAULT) -> dict:
    betas = linear_beta_schedule(T)
    return precompute_schedule_tensors(betas)


def get_m1_schedule() -> dict:
    betas = linear_beta_schedule(N_DIFFUSION_STEPS_M1)
    return precompute_schedule_tensors(betas)
