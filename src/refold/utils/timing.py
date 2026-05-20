"""Timing utilities."""

from __future__ import annotations

import time
from contextlib import contextmanager
from typing import Generator


@contextmanager
def timer(label: str = "") -> Generator[None, None, None]:
    """Context manager that prints elapsed time."""
    t0 = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - t0
        print(f"{label}: {elapsed:.2f}s" if label else f"Elapsed: {elapsed:.2f}s")


class Stopwatch:
    def __init__(self) -> None:
        self._start = time.perf_counter()

    def elapsed(self) -> float:
        return time.perf_counter() - self._start

    def reset(self) -> float:
        elapsed = self.elapsed()
        self._start = time.perf_counter()
        return elapsed
