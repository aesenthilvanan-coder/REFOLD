"""Load PhosphoSitePlus post-translational modification data."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from refold.constants import RAW_DIR, PROCESSED_DIR

logger = logging.getLogger(__name__)

PHOSPHOSITE_RAW_PATH = RAW_DIR / "phosphosite" / "phosphosite_plus.csv"
PHOSPHOSITE_PROCESSED_PATH = PROCESSED_DIR / "phosphosite_ptms.parquet"


def load_phosphosite_data(force: bool = False) -> pd.DataFrame:
    """
    Load PhosphoSitePlus PTM data.

    PhosphoSitePlus requires manual download from https://www.phosphosite.org/
    under a research license. Place the CSV at data/raw/phosphosite/phosphosite_plus.csv.
    """
    if PHOSPHOSITE_PROCESSED_PATH.exists() and not force:
        logger.info(f"Loading cached PhosphoSite data from {PHOSPHOSITE_PROCESSED_PATH}")
        return pd.read_parquet(PHOSPHOSITE_PROCESSED_PATH)

    if not PHOSPHOSITE_RAW_PATH.exists():
        logger.warning(
            f"PhosphoSitePlus data not found at {PHOSPHOSITE_RAW_PATH}. "
            "Download manually from https://www.phosphosite.org/staticDownloads"
        )
        return pd.DataFrame(columns=["uniprot_id", "position", "ptm_type", "modification"])

    logger.info("Parsing PhosphoSitePlus data...")
    df = _parse_phosphosite(PHOSPHOSITE_RAW_PATH)

    PHOSPHOSITE_PROCESSED_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(PHOSPHOSITE_PROCESSED_PATH, index=False)
    logger.info(f"PhosphoSite: {len(df):,} PTMs saved")
    return df


def _parse_phosphosite(path: Path) -> pd.DataFrame:
    """Parse PhosphoSitePlus CSV."""
    try:
        df = pd.read_csv(path, low_memory=False, on_bad_lines="skip", comment="#")
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

        rows = []
        for _, row in df.iterrows():
            uniprot = str(row.get("acc_id", row.get("uniprot_id", ""))).strip()
            mod_rsd = str(row.get("mod_rsd", "")).strip()
            ptm_type = str(row.get("modification", row.get("modification_type", ""))).strip()

            if not uniprot or not mod_rsd:
                continue

            # Parse position from mod_rsd (e.g., "S142-p" or "T123")
            try:
                pos_str = "".join(c for c in mod_rsd if c.isdigit())
                pos = int(pos_str) if pos_str else None
            except ValueError:
                pos = None

            if pos is None:
                continue

            rows.append({
                "uniprot_id": uniprot,
                "position": pos,
                "modification": mod_rsd,
                "ptm_type": ptm_type,
            })

        return pd.DataFrame(rows) if rows else pd.DataFrame(
            columns=["uniprot_id", "position", "modification", "ptm_type"]
        )
    except Exception as e:
        logger.error(f"PhosphoSite parsing failed: {e}")
        return pd.DataFrame(columns=["uniprot_id", "position", "modification", "ptm_type"])


def get_ptms_for_protein(uniprot_id: str, df: pd.DataFrame) -> list[dict]:
    """Return list of PTMs for a given UniProt ID."""
    sub = df[df["uniprot_id"] == uniprot_id]
    return sub.to_dict("records")
