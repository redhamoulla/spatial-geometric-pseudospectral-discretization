PYTHON ?= python3

.PHONY: test reproduce verify
test:
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v

reproduce:
	$(PYTHON) experiments/reproduce.py --check

verify:
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v
	$(PYTHON) experiments/reproduce.py --check --no-figures --output /tmp/gpps-verify
