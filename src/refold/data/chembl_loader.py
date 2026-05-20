"""Load and preprocess ChEMBL binding affinity data."""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Optional

import pandas as pd

from refold.constants import RAW_DIR, PROCESSED_DIR, CHEMBL_URL

logger = logging.getLogger(__name__)

CHEMBL_RAW_PATH = RAW_DIR / "chembl" / "chembl.db"
CHEMBL_PROCESSED_PATH = PROCESSED_DIR / "chembl_binding.parquet"


def load_chembl_binding_data(
    force: bool = False,
    min_confidence: int = 5,
) -> pd.DataFrame:
    """
    Load ChEMBL compound-target binding data.

    Returns DataFrame with columns: smiles, uniprot_id, pchembl_value, activity_type.
    """
    if CHEMBL_PROCESSED_PATH.exists() and not force:
        logger.info(f"Loading cached ChEMBL data from {CHEMBL_PROCESSED_PATH}")
        return pd.read_parquet(CHEMBL_PROCESSED_PATH)

    if not CHEMBL_RAW_PATH.exists():
        logger.warning(
            f"ChEMBL database not found at {CHEMBL_RAW_PATH}. "
            f"Download with: wget {CHEMBL_URL}"
        )
        return pd.DataFrame(columns=["smiles", "uniprot_id", "pchembl_value", "activity_type"])

    logger.info("Querying ChEMBL database...")
    df = _query_chembl_db(CHEMBL_RAW_PATH, min_confidence)

    CHEMBL_PROCESSED_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(CHEMBL_PROCESSED_PATH, index=False)
    logger.info(f"ChEMBL: {len(df):,} binding measurements saved")
    return df


def _query_chembl_db(db_path: Path, min_confidence: int) -> pd.DataFrame:
    """Query relevant binding data from ChEMBL SQLite."""
    query = """
    SELECT DISTINCT
        cs.canonical_smiles AS smiles,
        td.chembl_id AS target_chembl_id,
        act.pchembl_value,
        act.standard_type AS activity_type,
        tc.accession AS uniprot_id
    FROM activities act
    JOIN assays ass ON act.assay_id = ass.assay_id
    JOIN target_dictionary td ON ass.tid = td.tid
    JOIN target_components tc2 ON td.tid = tc2.tid
    JOIN component_sequences tc ON tc2.component_id = tc.component_id
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    JOIN compound_structures cs ON md.molregno = cs.molregno
    WHERE
        ass.confidence_score >= ?
        AND act.pchembl_value IS NOT NULL
        AND tc.accession IS NOT NULL
        AND cs.canonical_smiles IS NOT NULL
        AND td.target_type = 'SINGLE PROTEIN'
        AND act.standard_type IN ('IC50', 'Ki', 'Kd', 'EC50')
    """
    try:
        conn = sqlite3.connect(str(db_path))
        df = pd.read_sql_query(query, conn, params=(min_confidence,))
        conn.close()
        return df
    except Exception as e:
        logger.error(f"ChEMBL query failed: {e}")
        return pd.DataFrame(columns=["smiles", "uniprot_id", "pchembl_value", "activity_type"])
