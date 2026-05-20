# REFOLD: Restorative Engine for Folding Optimization via Ligand Design

REFOLD is a fully automated, end-to-end computational pipeline that takes any human disease-associated missense mutation as input and outputs:

1. A binary prediction of whether the mutation causes misfolding (amenable to pharmacological chaperone rescue)
2. The 3D structure of the mutant protein in its partially unfolded state with transiently exposed binding pockets
3. De novo designed small molecule candidates predicted to bind those transient pockets and thermodynamically stabilize the correct fold
4. Predicted pharmacokinetic properties, synthesizability scores, and cellular rescue probability for each candidate

## Quick Start

```bash
# Install
pip install -e ".[ml]"

# Run on Fabry disease mutation (GLA Y152C)
refold-single \
    --uniprot P06280 \
    --position 152 \
    --wildtype Y \
    --mutant C \
    --gene GLA \
    --disease "Fabry disease"
```

## Installation

```bash
# Core dependencies
pip install -e .

# ML dependencies (ESM-2, torch-geometric)
pip install -e ".[ml]"
pip install fair-esm torch-geometric

# External tools (macOS)
brew install fpocket openbabel dssp
```

## Data Setup

```bash
# Download raw data (ClinVar, ThermoMutDB)
make download

# Preprocess ClinVar (~30-60 min)
make preprocess
```

## Training

```bash
# Train rescue classifier
make train-rescue

# Train molecule generator (uses synthetic data by default)
make train-generator
```

## Pipeline Stages

| Stage | Description | Time/mutation |
|-------|-------------|---------------|
| 1 | Classify only (ΔΔG + rescue probability) | ~10s |
| 2 | Stage 1 + pocket detection (ANM + fpocket) | ~2min |
| 3 | Stage 2 + molecule generation | ~10min |

## Architecture

- **ESM-1v**: ΔΔG prediction via masked marginal probability
- **ESM-2 (150M)**: Protein embeddings for rescue classifier input
- **ANM/ENM**: Conformational sampling (50 states per mutation)
- **fpocket**: Alpha-sphere binding pocket detection
- **DDPM Diffusion**: Pocket-conditioned de novo 3D molecule generation
- **RescueClassifier**: GNN + ESM-2 + thermodynamic + evolutionary features

## Benchmark Targets

| Metric | Target |
|--------|--------|
| AUROC  | ≥ 0.85 |
| AUPRC  | ≥ 0.70 |
| MCC    | ≥ 0.55 |
| F1     | ≥ 0.65 |
| ΔΔG Pearson r | ≥ 0.60 |

```bash
# Run benchmark evaluation
make benchmark
```

## Known Chaperone-Responsive Benchmark Variants

- **Fabry disease** (GLA): Y152C, R215W, R231H
- **Pompe disease** (GAA): L525S, P600L
- **Gaucher disease** (GBA): N444S, L296V
- **Cystic fibrosis** (CFTR): R551H
- **Transthyretin amyloidosis** (TTR): V30M, V122I
- **Li-Fraumeni syndrome** (TP53): R175H, R248W, R249S

## License

MIT
