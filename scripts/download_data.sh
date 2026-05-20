#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
RAW_DIR="$ROOT_DIR/data/raw"

echo "Creating data directories..."
mkdir -p "$RAW_DIR"/{clinvar,alphafold,thermomutdb,protherm,chembl,bindingdb,phosphosite}

# ClinVar XML (~2.5GB)
echo "Downloading ClinVar XML..."
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
CLINVAR_OUT="$RAW_DIR/clinvar/ClinVarFullRelease_00-latest.xml.gz"
if [ ! -f "$CLINVAR_OUT" ]; then
    wget -c -O "$CLINVAR_OUT" "$CLINVAR_URL" && echo "ClinVar downloaded." || echo "ClinVar download failed."
else
    echo "ClinVar already exists."
fi

# UniProt human proteome FASTA
echo "Downloading UniProt human FASTA..."
UNIPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
if [ ! -f "$RAW_DIR/uniprot_human.fasta.gz" ]; then
    wget -c -O "$RAW_DIR/uniprot_human.fasta.gz" "$UNIPROT_URL" && echo "UniProt FASTA downloaded." || echo "UniProt download failed."
fi

# ThermoMutDB
echo "Downloading ThermoMutDB..."
THERMO_URL="https://biosig.lab.uq.edu.au/thermomutdb/static/media/thermomutdb_json.zip"
if [ ! -f "$RAW_DIR/thermomutdb/thermomutdb_json.zip" ]; then
    wget -c -O "$RAW_DIR/thermomutdb/thermomutdb_json.zip" "$THERMO_URL" && echo "ThermoMutDB downloaded." || echo "ThermoMutDB download failed."
fi

# AlphaFold bulk download (commented out — fetched on-demand during pipeline)
# echo "Downloading AlphaFold human proteome..."
# AF_URL="https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar"
# wget -c -O "$RAW_DIR/alphafold/UP000005640_9606_HUMAN_v4.tar" "$AF_URL"

echo "Data download complete."
