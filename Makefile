PYTHON ?= python3

.PHONY: test generate screen train-properties optimize report validate-manifests pipeline

test:
	$(PYTHON) -m unittest discover -s tests -q

generate:
	$(PYTHON) -m peptide.run_generation

screen:
	$(PYTHON) -m scoring.run_screening

train-properties:
	$(PYTHON) -m models.run_property_training

optimize:
	$(PYTHON) models/run_optimization.py --top-k 25

report:
	$(PYTHON) reports/run_report.py

validate-manifests:
	$(PYTHON) -m core.run_manifest_validation

pipeline: generate screen train-properties optimize report validate-manifests
