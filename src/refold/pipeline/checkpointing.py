"""Pipeline checkpointing utilities."""

from __future__ import annotations

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class PipelineCheckpoint:
    """Manages resumable pipeline progress."""

    def __init__(self, checkpoint_file: Path):
        self.checkpoint_file = checkpoint_file
        self._data: dict = {"completed": [], "failed": [], "metadata": {}}
        if checkpoint_file.exists():
            self._load()

    def _load(self) -> None:
        try:
            with open(self.checkpoint_file) as f:
                self._data = json.load(f)
        except Exception as e:
            logger.warning(f"Could not load checkpoint {self.checkpoint_file}: {e}")

    def _save(self) -> None:
        self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.checkpoint_file, "w") as f:
            json.dump(self._data, f, indent=2)

    def is_completed(self, key: str) -> bool:
        return key in self._data["completed"]

    def mark_completed(self, key: str) -> None:
        if key not in self._data["completed"]:
            self._data["completed"].append(key)
        self._save()

    def mark_failed(self, key: str, error: str = "") -> None:
        entry = {"key": key, "error": error}
        self._data["failed"].append(entry)
        self._save()

    @property
    def n_completed(self) -> int:
        return len(self._data["completed"])

    @property
    def n_failed(self) -> int:
        return len(self._data["failed"])
