"""Load and preprocess ProThermDB stability data."""

from __future__ import annotations

import logging
import zipfile
from io import BytesIO
from pathlib import Path

import pandas as pd
import requests

from refold.constants import RAW_DIR, PROCESSED_DIR, PROTHERM_URL

logger = logging.getLogger(__name__)

PROTHERM_RAW_PATH = RAW_DIR / "protherm" / "prothermdb_raw.csv"
PROTHERM_PROCESSED_PATH = PROCESSED_DIR / "protherm_stability.parquet"


def load_protherm_data(force: bool = False) -> pd.DataFrame:
    """Load ProThermDB stability measurements."""
    if PROTHERM_PROCESSED_PATH.exists() and not force:
        logger.info(f"Loading cached ProTherm data from {PROTHERM_PROCESSED_PATH}")
        return pd.read_parquet(PROTHERM_PROCESSED_PATH)

    if not PROTHERM_RAW_PATH.exists():
        logger.info("Downloading ProThermDB...")
        _download_protherm()

    if not PROTHERM_RAW_PATH.exists():
        logger.warning("ProThermDB not available")
        return pd.DataFrame(columns=["uniprot_id", "mutation", "ddg", "ph", "temperature"])

    logger.info("Parsing ProThermDB...")
    df = _parse_protherm(PROTHERM_RAW_PATH)

    PROTHERM_PROCESSED_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(PROTHERM_PROCESSED_PATH, index=False)
    logger.info(f"ProThermDB: {len(df):,} measurements saved")
    return df


def _download_protherm() -> None:
    PROTHERM_RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    try:
        response = requests.get(PROTHERM_URL, stream=True, timeout=120)
        response.raise_for_status()
        content = BytesIO(response.content)
        if PROTHERM_URL.endswith(".zip"):
            with zipfile.ZipFile(content) as zf:
                csv_name = next(
                    (n for n in zf.namelist() if n.endswith(".csv") or n.endswith(".tsv")),
                    zf.namelist()[0],
                )
                with zf.open(csv_name) as src, open(PROTHERM_RAW_PATH, "wb") as dst:
                    dst.write(src.read())
        else:
            with open(PROTHERM_RAW_PATH, "wb") as f:
                f.write(response.content)
        logger.info(f"ProThermDB downloaded to {PROTHERM_RAW_PATH}")
    except Exception as e:
        logger.error(f"ProThermDB download failed: {e}")


def _parse_protherm(path: Path) -> pd.DataFrame:
    """Parse ProThermDB CSV into normalized DataFrame."""
    try:
        # Try tab-separated first
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False, on_bad_lines="skip")
        except Exception:
            df = pd.read_csv(path, low_memory=False, on_bad_lines="skip")

        # Normalize column names
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

        result_rows = []
        for _, row in df.iterrows():
            ddg = _extract_ddg(row)
            if ddg is None or abs(ddg) > 20:
                continue

            ph = _safe_float(row.get("ph", 7.0))
            temp = _safe_float(row.get("temperature", 25.0))
            if ph is None or not (5.0 <= ph <= 9.0):
                continue
            if temp is None or not (15.0 <= temp <= 45.0):
                continue

            uniprot = str(row.get("uniprot", row.get("uniprot_id", ""))).strip()
            mutation = str(row.get("mutation", "")).strip()

            if not uniprot or not mutation:
                continue

            result_rows.append({
                "uniprot_id": uniprot,
                "mutation": mutation,
                "ddg": ddg,
                "ph": ph,
                "temperature": temp,
            })

        return pd.DataFrame(result_rows) if result_rows else pd.DataFrame(
            columns=["uniprot_id", "mutation", "ddg", "ph", "temperature"]
        )

    except Exception as e:
        logger.error(f"ProTherm parsing failed: {e}")
        return pd.DataFrame(columns=["uniprot_id", "mutation", "ddg", "ph", "temperature"])


def _extract_ddg(row) -> float | None:
    for col in ["ddg", "delta_delta_g", "deltadeltag", "ddg_kcal_mol"]:
        val = row.get(col)
        if pd.notna(val):
            try:
                return float(val)
            except (ValueError, TypeError):
                pass
    return None


def _safe_float(val) -> float | None:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None
