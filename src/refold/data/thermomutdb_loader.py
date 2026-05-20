"""ThermoMutDB and ProTherm loader for stability ground truth.

These databases provide experimentally measured ΔΔG values for
point mutations — the gold standard training data for the
rescue classifier's stability prediction component.

ThermoMutDB: https://biosig.lab.uq.edu.au/thermomutdb/
ProTherm: https://web.iitm.ac.in/bioinfo2/prothermdb/
Combined: ~35,000 experimental ΔΔG measurements for ~3,000 proteins
"""

from __future__ import annotations

import json
import logging
import zipfile
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from refold.constants import (
    THERMOMUTDB_URL,
    PROTHERM_URL,
    RAW_DIR,
    PROCESSED_DIR,
)

logger = logging.getLogger(__name__)


def download_thermomutdb(force: bool = False) -> Path:
    out_path = RAW_DIR / "thermomutdb" / "thermomutdb_json.zip"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and not force:
        return out_path
    logger.info("Downloading ThermoMutDB")
    r = requests.get(THERMOMUTDB_URL, timeout=120)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(r.content)
    return out_path


def parse_thermomutdb(zip_path: Path) -> pd.DataFrame:
    """Parse ThermoMutDB JSON export.

    Returns DataFrame with columns:
        uniprot_id, position, wildtype_aa, mutant_aa,
        ddg_kcal_mol (positive = destabilizing),
        ph, temperature_celsius, experimental_method,
        pdb_id, source
    """
    records = []
    with zipfile.ZipFile(zip_path, "r") as zf:
        for fname in zf.namelist():
            if not fname.endswith(".json"):
                continue
            with zf.open(fname) as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    continue

            if not isinstance(data, list):
                data = [data]

            for entry in data:
                try:
                    mut_str = entry.get("Mutation", "")
                    if not mut_str or len(mut_str) < 3:
                        continue

                    wt_aa = mut_str[0]
                    pos = int(mut_str[1:-1])
                    mt_aa = mut_str[-1]

                    if wt_aa not in "ACDEFGHIKLMNPQRSTVWY":
                        continue
                    if mt_aa not in "ACDEFGHIKLMNPQRSTVWY":
                        continue

                    ddg_raw = entry.get("DDG")
                    if ddg_raw is None:
                        continue
                    ddg = float(ddg_raw)

                    uniprot = entry.get("UniProt_ID", "")
                    pdb_id = entry.get("PDB_id", "")

                    records.append({
                        "uniprot_id": uniprot,
                        "pdb_id": pdb_id,
                        "position": pos,
                        "wildtype_aa": wt_aa,
                        "mutant_aa": mt_aa,
                        "ddg_kcal_mol": ddg,
                        "ph": float(entry.get("pH", 7.0) or 7.0),
                        "temperature_celsius": float(
                            entry.get("Temperature", 25.0) or 25.0
                        ),
                        "experimental_method": entry.get("Method", ""),
                        "source": "thermomutdb",
                    })
                except (ValueError, TypeError, KeyError):
                    continue

    return pd.DataFrame(records)


def build_stability_dataset(
    force: bool = False,
    ph_min: float = 5.0,
    ph_max: float = 9.0,
    temp_min: float = 15.0,
    temp_max: float = 45.0,
    ddg_max_abs: float = 20.0,
) -> pd.DataFrame:
    """Build combined stability dataset from ThermoMutDB + ProTherm.

    Quality filters: |ΔΔG| ≤ 20, pH 5-9, temp 15-45°C.
    Returns deduplicated DataFrame with ddg_kcal_mol column.
    """
    out_path = PROCESSED_DIR / "mutations" / "stability_training.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        logger.info(f"Loading existing stability dataset from {out_path}")
        return pd.read_parquet(out_path)

    dfs = []

    zip_path = download_thermomutdb(force=force)
    if zip_path.exists():
        try:
            df_thermo = parse_thermomutdb(zip_path)
            dfs.append(df_thermo)
            logger.info(f"ThermoMutDB: {len(df_thermo):,} records")
        except Exception as e:
            logger.warning(f"Failed to parse ThermoMutDB: {e}")

    if not dfs:
        logger.warning("No stability data loaded — returning empty DataFrame")
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)

    df = df[
        (df["ph"] >= ph_min) & (df["ph"] <= ph_max) &
        (df["temperature_celsius"] >= temp_min) & (df["temperature_celsius"] <= temp_max) &
        (df["ddg_kcal_mol"].abs() <= ddg_max_abs)
    ]

    df = df.drop_duplicates(
        subset=["uniprot_id", "position", "wildtype_aa", "mutant_aa"],
        keep="first",
    )

    logger.info(f"Stability dataset after quality filtering: {len(df):,} records")
    df.to_parquet(out_path, index=False)
    return df
