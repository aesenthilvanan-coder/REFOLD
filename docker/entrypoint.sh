#!/bin/bash
set -e

echo "REFOLD container starting..."
echo "Device: ${REFOLD_DEVICE:-auto}"

# Create directory structure
mkdir -p data/raw/{clinvar,alphafold,thermomutdb,protherm,chembl,bindingdb}
mkdir -p data/processed/{mutations,structures,pockets,molecules,splits,id_mapping}
mkdir -p data/results/{rescue_predictions,pocket_predictions,molecule_candidates,clinvar_scan}
mkdir -p checkpoints/{rescue_classifier,pocket_detector,molecule_generator,admet_predictor}
mkdir -p logs/{training,inference}

exec "$@"
