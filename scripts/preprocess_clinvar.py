#!/usr/bin/env python3
"""Preprocess ClinVar XML to extract pathogenic missense mutations."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.data.clinvar_parser import download_clinvar_xml, parse_clinvar_to_mutations
from refold.data.uniprot_mapper import UniProtMapper
from refold.constants import RAW_DIR, PROCESSED_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess ClinVar XML")
    parser.add_argument("--clinvar-xml", type=Path, default=None,
                        help="Path to ClinVar XML.gz (downloads if absent)")
    parser.add_argument("--output", type=Path,
                        default=PROCESSED_DIR / "mutations" / "clinvar_missense.parquet")
    parser.add_argument("--force", action="store_true", help="Re-download and reprocess")
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("Step 1: Building UniProt ID mapper...")
    mapper = UniProtMapper()

    logger.info("Step 2: Downloading ClinVar XML (if needed)...")
    xml_path = args.clinvar_xml or download_clinvar_xml(force=args.force)

    logger.info("Step 3: Parsing ClinVar XML...")
    df = parse_clinvar_to_mutations(
        xml_path=xml_path,
        uniprot_mapper=mapper,
        output_path=args.output,
    )

    logger.info(f"Done. {len(df):,} pathogenic missense mutations saved to {args.output}")
    print(df.head())
    print(f"\nTotal: {len(df):,} mutations")
    print(f"Unique proteins: {df['uniprot_id'].nunique():,}")


if __name__ == "__main__":
    main()
