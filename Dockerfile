# Fully specified environment for maximum reproducibility (see docs/reproducibility.md).
# Build:  docker build -t outbreak-simulator .
# Run:    docker run --rm outbreak-simulator python examples/01_run_choir_outbreak.py
# Test:   docker run --rm outbreak-simulator pytest tests/ -q

FROM python:3.12.3-slim

WORKDIR /app

# System dependencies for matplotlib (headless rendering) and building any
# packages with native extensions.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt requirements-dev.txt pyproject.toml README.md ./
COPY src/ ./src/

RUN pip install --no-cache-dir -r requirements-dev.txt \
    && pip install --no-cache-dir -e .

COPY tests/ ./tests/
COPY examples/ ./examples/
COPY docs/ ./docs/
COPY scripts/ ./scripts/

# Headless matplotlib backend by default (no display in a container)
ENV MPLBACKEND=Agg

# Default: run the test suite, so `docker run outbreak-simulator` alone
# gives immediate confirmation the environment is correctly reproduced.
CMD ["pytest", "tests/", "-q"]
