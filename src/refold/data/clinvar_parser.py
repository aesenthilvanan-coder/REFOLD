"""ClinVar XML parser for extracting pathogenic missense mutations.

Downloads and parses the full ClinVar XML release.
Extracts all pathogenic/likely pathogenic missense variants
mapping to canonical human UniProt accessions.

Zero filtering based on disease — we want everything.
Output: parquet file at data/processed/mutations/clinvar_missense.parquet
"""

from __future__ import annotations

import gzip
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterator

import pandas as pd
import requests
from tqdm import tqdm

from refold.constants import (
    CLINVAR_XML_URL,
    RAW_DIR,
    PROCESSED_DIR,
    STANDARD_AAS,
)
from refold.types import Mutation
from refold.data.uniprot_mapper import UniProtMapper

logger = logging.getLogger(__name__)

PATHOGENIC_SIGNIFICANCES: frozenset[str] = frozenset({
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Pathogenic, low penetrance",
    "Likely pathogenic, low penetrance",
})


def download_clinvar_xml(force: bool = False) -> Path:
    """Download the full ClinVar XML release. File is ~2.5GB compressed. Resume-capable."""
    out_path = RAW_DIR / "clinvar" / "ClinVarFullRelease_00-latest.xml.gz"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        logger.info(f"ClinVar XML already exists at {out_path}")
        return out_path

    logger.info(f"Downloading ClinVar XML from {CLINVAR_XML_URL}")
    with requests.get(CLINVAR_XML_URL, stream=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(out_path, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc="ClinVar XML"
        ) as pbar:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
                pbar.update(len(chunk))

    logger.info(f"Downloaded ClinVar XML to {out_path}")
    return out_path


def _parse_hgvs_protein(hgvs_str: str) -> tuple[str, int, str] | None:
    """Parse a protein HGVS string of the form p.Xxx123Yyy.
    Returns (wildtype_aa_1letter, position, mutant_aa_1letter) or None.

    Examples:
        "p.Gly12Val" -> ("G", 12, "V")
        "p.Arg506Gln" -> ("R", 506, "Q")
    """
    from Bio.SeqUtils import seq1

    match = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvs_str)
    if not match:
        return None

    wt_3 = match.group(1)
    pos = int(match.group(2))
    mt_3 = match.group(3)

    if wt_3 == "Ter" or mt_3 == "Ter":
        return None

    try:
        wt_1 = seq1(wt_3)
        mt_1 = seq1(mt_3)
    except Exception:
        return None

    if wt_1 not in STANDARD_AAS or mt_1 not in STANDARD_AAS:
        return None

    if wt_1 == mt_1:
        return None

    return wt_1, pos, mt_1


def _iter_clinvar_records(xml_path: Path) -> Iterator[dict]:
    """Streaming XML parser for ClinVar records.
    Yields dicts with keys: variation_id, gene_symbol, clinical_significance,
    protein_change (HGVS), condition_name, rsid
    """
    opener = gzip.open if str(xml_path).endswith(".gz") else open

    with opener(xml_path, "rb") as f:
        context = ET.iterparse(f, events=("start", "end"))
        current_record: dict = {}
        inside_clinvar_set = False

        for event, elem in context:
            if event == "start" and elem.tag == "ClinVarSet":
                current_record = {}
                inside_clinvar_set = True

            if not inside_clinvar_set:
                continue

            if event == "end" and elem.tag == "ClinVarSet":
                yield current_record
                elem.clear()
                inside_clinvar_set = False
                current_record = {}

            elif event == "end":
                tag = elem.tag

                if tag == "ClinVarAssertion":
                    pass

                elif tag == "MeasureSet":
                    var_id = elem.get("ID")
                    if var_id:
                        current_record["variation_id"] = var_id

                elif tag == "Symbol":
                    elem_type = elem.find("ElementValue")
                    if elem_type is not None and elem_type.get("Type") == "Preferred":
                        current_record["gene_symbol"] = elem_type.text or ""

                elif tag == "AttributeSet":
                    attr = elem.find("Attribute")
                    if attr is not None:
                        attr_type = attr.get("Type", "")
                        if "protein change" in attr_type.lower():
                            current_record["protein_change"] = attr.text or ""

                elif tag == "ClinicalSignificance":
                    desc = elem.find("Description")
                    if desc is not None:
                        current_record["clinical_significance"] = desc.text or ""

                elif tag == "TraitSet":
                    trait_name = elem.find(".//*[@Type='Preferred']")
                    if trait_name is not None and "condition_name" not in current_record:
                        current_record["condition_name"] = trait_name.text or ""

                elif tag == "XRef":
                    if elem.get("DB") == "dbSNP":
                        current_record["rsid"] = f"rs{elem.get('ID', '')}"


def parse_clinvar_to_mutations(
    xml_path: Path,
    uniprot_mapper: UniProtMapper,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Parse ClinVar XML and extract all pathogenic missense mutations
    with UniProt accession mapping.

    This is the primary data ingestion function.
    Takes ~30-60 minutes for the full ClinVar release.

    Returns DataFrame with columns:
        uniprot_id, gene_symbol, position, wildtype_aa, mutant_aa,
        hgvs, clinvar_variation_id, clinical_significance,
        disease_name, rsid
    """
    if output_path is None:
        output_path = PROCESSED_DIR / "mutations" / "clinvar_missense.parquet"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        logger.info(f"Loading existing ClinVar data from {output_path}")
        return pd.read_parquet(output_path)

    logger.info("Parsing ClinVar XML — this will take 30–60 minutes")

    records = []
    n_total = 0
    n_pathogenic = 0
    n_missense = 0
    n_mapped = 0

    for record in tqdm(_iter_clinvar_records(xml_path), desc="Parsing ClinVar", unit="records"):
        n_total += 1

        sig = record.get("clinical_significance", "")
        if sig not in PATHOGENIC_SIGNIFICANCES:
            continue
        n_pathogenic += 1

        protein_change = record.get("protein_change", "")
        if not protein_change:
            continue

        parsed = _parse_hgvs_protein(protein_change)
        if parsed is None:
            continue
        wt_aa, position, mt_aa = parsed
        n_missense += 1

        gene_symbol = record.get("gene_symbol", "")
        uniprot_id = uniprot_mapper.gene_to_uniprot(gene_symbol)
        if uniprot_id is None:
            continue
        n_mapped += 1

        records.append({
            "uniprot_id": uniprot_id,
            "gene_symbol": gene_symbol,
            "position": position,
            "wildtype_aa": wt_aa,
            "mutant_aa": mt_aa,
            "hgvs": f"p.{wt_aa}{position}{mt_aa}",
            "clinvar_variation_id": record.get("variation_id"),
            "clinical_significance": sig,
            "disease_name": record.get("condition_name", ""),
            "rsid": record.get("rsid", ""),
        })

    logger.info(
        f"ClinVar parsing complete:\n"
        f"  Total records: {n_total:,}\n"
        f"  Pathogenic: {n_pathogenic:,}\n"
        f"  Missense: {n_missense:,}\n"
        f"  Mapped to UniProt: {n_mapped:,}"
    )

    df = pd.DataFrame(records)
    df = df.drop_duplicates(
        subset=["uniprot_id", "position", "wildtype_aa", "mutant_aa"],
        keep="first",
    )

    logger.info(f"Final deduplicated mutation count: {len(df):,}")
    df.to_parquet(output_path, index=False)
    logger.info(f"Saved to {output_path}")

    return df


def load_mutations_as_objects(df: pd.DataFrame) -> list[Mutation]:
    """Convert DataFrame rows to Mutation dataclass instances."""
    mutations = []
    for _, row in df.iterrows():
        try:
            m = Mutation(
                uniprot_id=row["uniprot_id"],
                position=int(row["position"]),
                wildtype_aa=row["wildtype_aa"],
                mutant_aa=row["mutant_aa"],
                gene_name=row.get("gene_symbol", ""),
                clinvar_id=str(row.get("clinvar_variation_id", "")),
                disease=row.get("disease_name", ""),
                source="clinvar",
            )
            mutations.append(m)
        except (AssertionError, ValueError) as e:
            logger.warning(f"Skipping malformed mutation row: {e}")
    return mutations
