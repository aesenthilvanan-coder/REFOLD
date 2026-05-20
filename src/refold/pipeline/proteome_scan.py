"""Proteome-scale scan across all human proteins."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from refold.constants import RESULTS_DIR

logger = logging.getLogger(__name__)


def run_proteome_scan(
    uniprot_ids: list[str] | None = None,
    stage: int = 1,
    n_workers: int = 4,
    device: "torch.device | None" = None,
    output_dir: Path | None = None,
) -> None:
    """Run REFOLD scan across a list of UniProt IDs and all their ClinVar variants."""
    from refold.utils.logging import setup_logging
    setup_logging()

    if output_dir is None:
        output_dir = RESULTS_DIR / "proteome_scan"
    output_dir.mkdir(parents=True, exist_ok=True)

    if uniprot_ids is None:
        from refold.data.uniprot_mapper import UniProtMapper
        mapper = UniProtMapper()
        uniprot_ids = mapper.all_human_uniprot_ids
        logger.info(f"Scanning all {len(uniprot_ids):,} human UniProt entries")

    from refold.constants import PROCESSED_DIR
    mutations_path = PROCESSED_DIR / "mutations" / "clinvar_missense.parquet"

    if not mutations_path.exists():
        raise FileNotFoundError(
            "ClinVar mutations not found. Run: python scripts/preprocess_clinvar.py"
        )

    import pandas as pd
    df = pd.read_parquet(mutations_path)
    df = df[df["uniprot_id"].isin(set(uniprot_ids))]
    logger.info(f"Found {len(df):,} mutations for target proteins")

    scan_path = output_dir / "mutations_subset.parquet"
    df.to_parquet(scan_path, index=False)

    from refold.pipeline.batch_clinvar import run_clinvar_scan
    run_clinvar_scan(
        mutations_parquet=scan_path,
        stage=stage,
        device=device,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="REFOLD proteome scan")
    parser.add_argument("--uniprot-list", type=Path, default=None,
                        help="File with UniProt IDs, one per line")
    parser.add_argument("--stage", type=int, default=1, choices=[1, 2, 3])
    parser.add_argument("--device", default=None)
    args = parser.parse_args()

    import torch
    device = torch.device(args.device) if args.device else None

    uids = None
    if args.uniprot_list:
        with open(args.uniprot_list) as f:
            uids = [line.strip() for line in f if line.strip()]

    run_proteome_scan(uniprot_ids=uids, stage=args.stage, device=device)


if __name__ == "__main__":
    main()
