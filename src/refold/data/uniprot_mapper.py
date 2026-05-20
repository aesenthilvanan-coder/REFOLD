"""UniProt ID mapping utilities.

Maps between:
- Gene symbols → canonical UniProt accessions
- Entrez gene IDs → UniProt
- Ensembl gene IDs → UniProt
- HGNC IDs → UniProt

Uses the UniProt ID mapping API and local cache.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import requests

from refold.constants import PROCESSED_DIR

logger = logging.getLogger(__name__)

MAPPING_CACHE_PATH: Path = PROCESSED_DIR / "id_mapping" / "uniprot_mapping.parquet"


class UniProtMapper:
    """Bidirectional ID mapper between gene identifiers and UniProt accessions.
    Downloads and caches the full human proteome ID mapping on first use.
    """

    def __init__(self, cache_path: Path = MAPPING_CACHE_PATH):
        self.cache_path = cache_path
        self._gene_to_uniprot: dict[str, str] = {}
        self._uniprot_to_gene: dict[str, str] = {}
        self._entrez_to_uniprot: dict[str, str] = {}
        self._ensembl_to_uniprot: dict[str, str] = {}
        self._loaded = False

    def _load(self) -> None:
        if self._loaded:
            return
        if self.cache_path.exists():
            self._load_from_cache()
        else:
            self._download_and_build_cache()
        self._loaded = True
        logger.info(
            f"UniProt mapper loaded: {len(self._gene_to_uniprot):,} genes, "
            f"{len(self._uniprot_to_gene):,} UniProt accessions"
        )

    def _download_and_build_cache(self) -> None:
        """Download human proteome ID mapping from UniProt API.
        Builds mapping tables and saves to parquet cache.
        """
        logger.info("Downloading UniProt ID mapping for human proteome")
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)

        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
            "format": "tsv",
            "query": "organism_id:9606 AND reviewed:true",
            "fields": "accession,gene_names,gene_oln,xref_geneid,xref_ensembl",
        }

        r = requests.get(url, params=params, timeout=120)
        r.raise_for_status()

        lines = r.text.strip().split("\n")
        records = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            accession = parts[0].strip()
            gene_names_raw = parts[1].strip()
            geneid_raw = parts[3].strip()
            ensembl_raw = parts[4].strip()

            gene_names = []
            if gene_names_raw:
                for token in gene_names_raw.split():
                    gene_names.append(token.upper())

            entrez_ids = []
            if geneid_raw:
                for eid in geneid_raw.split(";"):
                    eid = eid.strip()
                    if eid:
                        entrez_ids.append(eid)

            ensembl_ids = []
            if ensembl_raw:
                for ensg in ensembl_raw.split(";"):
                    ensg = ensg.strip().split("[")[0].strip()
                    if ensg.startswith("ENSG"):
                        ensembl_ids.append(ensg)

            primary_gene = gene_names[0] if gene_names else ""
            records.append({
                "uniprot_id": accession,
                "gene_symbol": primary_gene,
                "all_gene_names": "|".join(gene_names),
                "entrez_ids": "|".join(entrez_ids),
                "ensembl_ids": "|".join(ensembl_ids),
            })

        df = pd.DataFrame(records)
        df.to_parquet(self.cache_path, index=False)
        self._build_indices(df)

    def _load_from_cache(self) -> None:
        df = pd.read_parquet(self.cache_path)
        self._build_indices(df)

    def _build_indices(self, df: pd.DataFrame) -> None:
        for _, row in df.iterrows():
            uid = row["uniprot_id"]
            for gene in str(row["all_gene_names"]).split("|"):
                if gene:
                    self._gene_to_uniprot[gene.upper()] = uid
            if row["gene_symbol"]:
                self._uniprot_to_gene[uid] = row["gene_symbol"]
            for eid in str(row["entrez_ids"]).split("|"):
                if eid:
                    self._entrez_to_uniprot[eid] = uid
            for ensg in str(row["ensembl_ids"]).split("|"):
                if ensg:
                    self._ensembl_to_uniprot[ensg] = uid

    def gene_to_uniprot(self, gene_symbol: str) -> str | None:
        """Map a gene symbol (any case) to canonical UniProt accession."""
        self._load()
        return self._gene_to_uniprot.get(gene_symbol.upper())

    def uniprot_to_gene(self, uniprot_id: str) -> str | None:
        """Map UniProt accession to primary gene symbol."""
        self._load()
        return self._uniprot_to_gene.get(uniprot_id)

    def entrez_to_uniprot(self, entrez_id: str) -> str | None:
        self._load()
        return self._entrez_to_uniprot.get(str(entrez_id))

    def ensembl_to_uniprot(self, ensembl_id: str) -> str | None:
        self._load()
        return self._ensembl_to_uniprot.get(ensembl_id)

    def is_valid_uniprot(self, uniprot_id: str) -> bool:
        self._load()
        return uniprot_id in self._uniprot_to_gene

    @property
    def all_human_uniprot_ids(self) -> list[str]:
        self._load()
        return list(self._uniprot_to_gene.keys())
