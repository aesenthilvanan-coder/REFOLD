"""AlphaFold structure fetcher and local cache manager.

Downloads AlphaFold v4 predicted structures for human proteins.
Maintains a local cache keyed by UniProt accession.
Supports batch download with resume capability.
"""

from __future__ import annotations

import gzip
import logging
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from tqdm import tqdm

from refold.constants import ALPHAFOLD_BASE_URL, RAW_DIR

logger = logging.getLogger(__name__)

ALPHAFOLD_VERSION: str = "v4"
AF_STRUCTURE_DIR: Path = RAW_DIR / "alphafold"
AF_STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


def _af_url(uniprot_id: str) -> str:
    """Construct AlphaFold download URL for a UniProt accession."""
    fname = f"AF-{uniprot_id}-F1-model_{ALPHAFOLD_VERSION}.pdb.gz"
    return f"{ALPHAFOLD_BASE_URL}{fname}"


def _af_local_path(uniprot_id: str, compressed: bool = False) -> Path:
    ext = ".pdb.gz" if compressed else ".pdb"
    return AF_STRUCTURE_DIR / f"AF-{uniprot_id}-F1-model_{ALPHAFOLD_VERSION}{ext}"


def fetch_structure(
    uniprot_id: str,
    force: bool = False,
    max_retries: int = 3,
    retry_delay: float = 2.0,
) -> Path | None:
    """Fetch and decompress an AlphaFold structure for a single UniProt ID.
    Returns the local .pdb path, or None if unavailable.

    Retry logic handles transient network failures.
    """
    local_pdb = _af_local_path(uniprot_id, compressed=False)

    if local_pdb.exists() and not force:
        return local_pdb

    url = _af_url(uniprot_id)
    compressed_path = _af_local_path(uniprot_id, compressed=True)

    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                logger.debug(f"No AlphaFold structure for {uniprot_id}")
                return None
            r.raise_for_status()
            with open(compressed_path, "wb") as f:
                f.write(r.content)
            break
        except requests.RequestException as e:
            if attempt == max_retries - 1:
                logger.warning(
                    f"Failed to fetch {uniprot_id} after {max_retries} attempts: {e}"
                )
                return None
            time.sleep(retry_delay * (attempt + 1))

    try:
        with gzip.open(compressed_path, "rb") as gz:
            with open(local_pdb, "wb") as pdb:
                pdb.write(gz.read())
        compressed_path.unlink()
        return local_pdb
    except Exception as e:
        logger.warning(f"Failed to decompress {compressed_path}: {e}")
        return None


def batch_fetch_structures(
    uniprot_ids: list[str],
    force: bool = False,
    n_workers: int = 8,
    rate_limit_per_second: float = 5.0,
) -> dict[str, Path | None]:
    """Fetch AlphaFold structures for a list of UniProt IDs.

    Rate-limited to respect EBI servers.
    Returns dict mapping uniprot_id → local Path (or None if unavailable).
    """
    results: dict[str, Path | None] = {}
    semaphore = threading.Semaphore(n_workers)
    rate_lock = threading.Lock()
    last_request_time = [0.0]
    min_interval = 1.0 / rate_limit_per_second

    def fetch_with_rate_limit(uid: str) -> tuple[str, Path | None]:
        with semaphore:
            with rate_lock:
                now = time.time()
                elapsed = now - last_request_time[0]
                if elapsed < min_interval:
                    time.sleep(min_interval - elapsed)
                last_request_time[0] = time.time()
            return uid, fetch_structure(uid, force=force)

    to_fetch = [
        uid for uid in uniprot_ids
        if force or not _af_local_path(uid).exists()
    ]
    for uid in uniprot_ids:
        if _af_local_path(uid).exists() and not force:
            results[uid] = _af_local_path(uid)

    if not to_fetch:
        logger.info("All structures already cached")
        return results

    logger.info(f"Fetching {len(to_fetch)} AlphaFold structures")
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(fetch_with_rate_limit, uid): uid for uid in to_fetch}
        for future in tqdm(as_completed(futures), total=len(to_fetch), desc="AlphaFold fetch"):
            uid, path = future.result()
            results[uid] = path

    n_success = sum(1 for v in results.values() if v is not None)
    logger.info(
        f"Fetched {n_success}/{len(uniprot_ids)} structures "
        f"({len(uniprot_ids) - n_success} unavailable)"
    )
    return results


def get_structure_path(uniprot_id: str) -> Path | None:
    """Get cached structure path without downloading."""
    p = _af_local_path(uniprot_id)
    return p if p.exists() else None


def structure_exists(uniprot_id: str) -> bool:
    """Check if structure is locally cached."""
    return _af_local_path(uniprot_id).exists()
