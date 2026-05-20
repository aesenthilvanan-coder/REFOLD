"""I/O utilities for REFOLD."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""

    def default(self, obj: Any) -> Any:
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def save_result_json(result_dict: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(result_dict, f, indent=2, cls=NumpyEncoder)
    logger.info(f"Saved result JSON to {path}")


def save_molecules_csv(molecules: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(molecules)
    df.to_csv(path, index=False)
    logger.info(f"Saved {len(molecules)} molecules to {path}")


def load_config(config_path: Path) -> dict[str, Any]:
    import yaml

    with open(config_path) as f:
        return yaml.safe_load(f)
