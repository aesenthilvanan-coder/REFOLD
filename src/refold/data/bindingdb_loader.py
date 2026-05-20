"""Load and preprocess BindingDB data."""

from __future__ import annotations

import logging
import zipfile
from io import BytesIO
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from refold.constants import RAW_DIR, PROCESSED_DIR, BINDINGDB_URL

logger = logging.getLogger(__name__)

BINDINGDB_RAW_PATH = RAW_DIR / "bindingdb" / "BindingDB_All.tsv"
BINDINGDB_PROCESSED_PATH = PROCESSED_DIR / "bindingdb_binding.parquet"

_COLUMNS_OF_INTEREST = [
    "Ligand SMILES",
    "Target Name",
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "UniProt (SwissProt) Primary ID of Target Chain",
]


def load_bindingdb_data(force: bool = False) -> pd.DataFrame:
    """Load BindingDB binding data as a DataFrame."""
    if BINDINGDB_PROCESSED_PATH.exists() and not force:
        logger.info(f"Loading cached BindingDB data from {BINDINGDB_PROCESSED_PATH}")
        return pd.read_parquet(BINDINGDB_PROCESSED_PATH)

    if not BINDINGDB_RAW_PATH.exists():
        logger.info("Downloading BindingDB...")
        _download_bindingdb()

    if not BINDINGDB_RAW_PATH.exists():
        logger.warning("BindingDB file not available")
        return pd.DataFrame(columns=["smiles", "uniprot_id", "affinity_nm", "activity_type"])

    logger.info("Parsing BindingDB TSV...")
    df = _parse_bindingdb(BINDINGDB_RAW_PATH)

    BINDINGDB_PROCESSED_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(BINDINGDB_PROCESSED_PATH, index=False)
    logger.info(f"BindingDB: {len(df):,} records saved")
    return df


def _download_bindingdb() -> None:
    BINDINGDB_RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    try:
        response = requests.get(BINDINGDB_URL, stream=True, timeout=300)
        response.raise_for_status()
        content = BytesIO(response.content)
        with zipfile.ZipFile(content) as zf:
            tsv_name = next(n for n in zf.namelist() if n.endswith(".tsv"))
            with zf.open(tsv_name) as src, open(BINDINGDB_RAW_PATH, "wb") as dst:
                dst.write(src.read())
        logger.info(f"BindingDB downloaded to {BINDINGDB_RAW_PATH}")
    except Exception as e:
        logger.error(f"BindingDB download failed: {e}")


def _parse_bindingdb(path: Path) -> pd.DataFrame:
    """Parse BindingDB TSV into normalized DataFrame."""
    try:
        df = pd.read_csv(
            path, sep="\t", usecols=_COLUMNS_OF_INTEREST, low_memory=False,
            on_bad_lines="skip",
        )
    except Exception as e:
        logger.error(f"BindingDB parsing failed: {e}")
        return pd.DataFrame()

    df = df.rename(columns={
        "Ligand SMILES": "smiles",
        "UniProt (SwissProt) Primary ID of Target Chain": "uniprot_id",
    })
    df = df.dropna(subset=["smiles", "uniprot_id"])

    # Melt affinity columns
    rows = []
    for _, row in df.iterrows():
        for col, atype in [("Ki (nM)", "Ki"), ("IC50 (nM)", "IC50"),
                            ("Kd (nM)", "Kd"), ("EC50 (nM)", "EC50")]:
            val = row.get(col)
            if pd.notna(val) and str(val).strip():
                try:
                    nm = float(str(val).replace(">", "").replace("<", "").strip())
                    rows.append({
                        "smiles": row["smiles"],
                        "uniprot_id": row["uniprot_id"],
                        "affinity_nm": nm,
                        "activity_type": atype,
                    })
                    break
                except ValueError:
                    pass

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["smiles", "uniprot_id", "affinity_nm", "activity_type"]
    )
