"""Abstract base class for REFOLD pipeline stages."""

from __future__ import annotations

import logging
import time
from abc import ABC, abstractmethod
from typing import Any, Optional

from refold.types import Mutation, REFOLDResult

logger = logging.getLogger(__name__)


class BasePipelineStage(ABC):
    """A single stage in the REFOLD pipeline."""

    name: str = "unnamed_stage"

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    @abstractmethod
    def run(self, mutation: Mutation, context: dict[str, Any]) -> dict[str, Any]:
        """Execute this stage and return updated context dict."""
        ...

    def __call__(self, mutation: Mutation, context: dict[str, Any]) -> dict[str, Any]:
        if not self.enabled:
            logger.debug(f"Stage {self.name} skipped (disabled)")
            return context
        t0 = time.perf_counter()
        try:
            result = self.run(mutation, context)
            elapsed = time.perf_counter() - t0
            logger.debug(f"Stage {self.name} completed in {elapsed:.2f}s")
            return result
        except Exception as e:
            logger.error(f"Stage {self.name} failed: {e}")
            context["stage_errors"] = context.get("stage_errors", [])
            context["stage_errors"].append({"stage": self.name, "error": str(e)})
            return context


class BasePipeline(ABC):
    """Base class for the full REFOLD pipeline."""

    def __init__(self, device: Optional[Any] = None):
        self._device = device
        self._stages: list[BasePipelineStage] = []

    @abstractmethod
    def run(self, mutation: Mutation) -> REFOLDResult:
        """Execute the full pipeline for a single mutation."""
        ...

    def add_stage(self, stage: BasePipelineStage) -> None:
        self._stages.append(stage)

    @property
    def device(self) -> Any:
        if self._device is None:
            from refold.utils.device import get_device
            self._device = get_device()
        return self._device
