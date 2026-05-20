"""Memory management helpers."""

from __future__ import annotations

import gc
import logging

logger = logging.getLogger(__name__)


def cleanup_mps(every_n: int = 50, step: int = 0) -> None:
    """Run MPS memory cleanup every `every_n` steps."""
    if step % every_n != 0:
        return
    try:
        import torch

        if torch.backends.mps.is_available():
            torch.mps.empty_cache()
    except Exception:
        pass
    gc.collect()


def log_memory_usage(tag: str = "") -> None:
    """Log current memory usage."""
    try:
        import torch

        if torch.backends.mps.is_available():
            alloc = torch.mps.current_allocated_memory() / 1e9
            driver = torch.mps.driver_allocated_memory() / 1e9
            logger.debug(f"[{tag}] MPS allocated={alloc:.2f}GB driver={driver:.2f}GB")
        elif torch.cuda.is_available():
            alloc = torch.cuda.memory_allocated() / 1e9
            reserved = torch.cuda.memory_reserved() / 1e9
            logger.debug(f"[{tag}] CUDA allocated={alloc:.2f}GB reserved={reserved:.2f}GB")
    except Exception:
        pass
