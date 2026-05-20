.PHONY: install install-ml test test-unit test-integration test-quick lint format \
        download preprocess train-rescue train-generator train-all \
        run-single run-clinvar-stage1 run-clinvar-stage2 run-clinvar-full \
        benchmark dev-test clean clean-data clean-checkpoints

PYTHON := python3
PIP := pip3
UNIPROT ?= P06280
POS ?= 152
WT ?= Y
MT ?= C
GENE ?= GLA
DISEASE ?= "Fabry disease"

# ── Installation ─────────────────────────────────────────────────────────────

install:
	$(PIP) install -e .

install-ml:
	$(PIP) install -e ".[ml]"
	$(PIP) install fair-esm torch-geometric

# ── Testing ──────────────────────────────────────────────────────────────────

test:
	$(PYTHON) -m pytest tests/ -v --tb=short

test-unit:
	$(PYTHON) -m pytest tests/unit/ -v --tb=short

test-integration:
	$(PYTHON) -m pytest tests/integration/ -v --tb=short -m integration

test-quick:
	$(PYTHON) -m pytest tests/unit/ -x -q --tb=short

# ── Code Quality ─────────────────────────────────────────────────────────────

lint:
	$(PYTHON) -m flake8 src/ tests/ --max-line-length=100
	$(PYTHON) -m mypy src/refold/ --ignore-missing-imports

format:
	$(PYTHON) -m black src/ tests/ scripts/ --line-length=100
	$(PYTHON) -m isort src/ tests/ scripts/

# ── Data Pipeline ─────────────────────────────────────────────────────────────

download:
	bash scripts/download_data.sh

preprocess: preprocess-clinvar preprocess-thermodb preprocess-structures

preprocess-clinvar:
	$(PYTHON) scripts/preprocess_clinvar.py

preprocess-thermodb:
	$(PYTHON) scripts/preprocess_thermodb.py

preprocess-structures:
	$(PYTHON) scripts/preprocess_structures.py

# ── Training ─────────────────────────────────────────────────────────────────

train-rescue:
	$(PYTHON) scripts/train_rescue_classifier.py

train-generator:
	$(PYTHON) scripts/train_molecule_generator.py --synthetic-data

train-pocket-detector:
	$(PYTHON) scripts/train_pocket_detector.py

train-admet:
	$(PYTHON) scripts/train_admet_predictor.py

train-all: train-rescue train-generator train-pocket-detector train-admet

# ── Running ──────────────────────────────────────────────────────────────────

run-single:
	$(PYTHON) scripts/run_single.py \
		--uniprot $(UNIPROT) \
		--position $(POS) \
		--wildtype $(WT) \
		--mutant $(MT) \
		--gene $(GENE) \
		--disease $(DISEASE)

run-clinvar-stage1:
	caffeinate -i $(PYTHON) scripts/run_clinvar_scan.py --stage 1

run-clinvar-stage2:
	caffeinate -i $(PYTHON) scripts/run_clinvar_scan.py --stage 2

run-clinvar-full:
	caffeinate -i $(PYTHON) scripts/run_clinvar_scan.py --stage 3

run-proteome:
	caffeinate -i $(PYTHON) scripts/run_proteome_scan.py --stage 1

# ── Evaluation ────────────────────────────────────────────────────────────────

benchmark:
	$(PYTHON) scripts/evaluate_all.py

# ── Development ──────────────────────────────────────────────────────────────

dev-test:
	$(PYTHON) -m pytest tests/unit/ -x -q --tb=short \
		&& $(PYTHON) scripts/run_single.py \
			--uniprot P06280 --position 152 --wildtype Y --mutant C \
			--gene GLA --disease "Fabry disease" \
			--n-conformations 3 --n-molecules 10

# ── Cleanup ──────────────────────────────────────────────────────────────────

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	rm -rf .pytest_cache .mypy_cache .coverage htmlcov/
	rm -rf *.egg-info dist/ build/

clean-data:
	@echo "This will delete all downloaded and processed data!"
	@read -p "Are you sure? [y/N] " confirm && [ "$$confirm" = "y" ] && rm -rf data/ || true

clean-checkpoints:
	@echo "This will delete all model checkpoints!"
	@read -p "Are you sure? [y/N] " confirm && [ "$$confirm" = "y" ] && rm -rf checkpoints/ || true
