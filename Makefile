.PHONY: install install-dev test test-cov lint format evidence-tables examples docker-build docker-test clean

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test:
	pytest tests/ -q

test-cov:
	pytest tests/ --cov=outbreak_simulator --cov-report=term-missing

lint:
	ruff check src/ tests/

format:
	ruff format src/ tests/

evidence-tables:
	python scripts/generate_evidence_tables.py

examples:
	python examples/01_run_choir_outbreak.py
	python examples/02_compare_interventions.py
	python examples/03_sensitivity_analysis.py
	python examples/04_transmission_reconstruction.py

docker-build:
	docker build -t outbreak-simulator .

docker-test:
	docker run --rm outbreak-simulator pytest tests/ -q

clean:
	find . -type d -name __pycache__ -exec rm -rf {} +
	rm -rf .pytest_cache .coverage htmlcov *.egg-info build dist
	rm -f examples/*.png
