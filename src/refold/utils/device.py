"""Device selection and memory management utilities."""

from __future__ import annotations

import gc
import logging
from typing import Any

logger = logging.getLogger(__name__)


def get_device() -> "torch.device":
    """Return best available device: MPS > CUDA > CPU."""
    import torch

    if torch.backends.mps.is_available() and torch.backends.mps.is_built():
        device = torch.device("mps")
        logger.info("Using MPS (Apple Silicon) device")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
        logger.info(f"Using CUDA device: {torch.cuda.get_device_name(0)}")
    else:
        device = torch.device("cpu")
        logger.info("Using CPU device")
    return device


def free_memory(device: "torch.device | None" = None) -> None:
    """Release unused memory on the given device."""
    import torch

    gc.collect()
    if device is None:
        device = get_device()
    if device.type == "mps":
        torch.mps.empty_cache()
    elif device.type == "cuda":
        torch.cuda.empty_cache()


def get_memory_stats(device: "torch.device | None" = None) -> dict[str, float]:
    """Return memory statistics in GB for the given device."""
    import torch

    if device is None:
        device = get_device()

    stats: dict[str, float] = {}
    if device.type == "mps":
        stats["allocated_gb"] = torch.mps.current_allocated_memory() / 1e9
        stats["driver_allocated_gb"] = torch.mps.driver_allocated_memory() / 1e9
    elif device.type == "cuda":
        mem = torch.cuda.memory_stats(device)
        stats["allocated_gb"] = mem.get("allocated_bytes.all.current", 0) / 1e9
        stats["reserved_gb"] = mem.get("reserved_bytes.all.current", 0) / 1e9
        info = torch.cuda.mem_get_info(device)
        stats["free_gb"] = info[0] / 1e9
        stats["total_gb"] = info[1] / 1e9
    else:
        stats["allocated_gb"] = 0.0
    return stats
