# Reproducibility Guide

## Reproducibility level claimed

Per the reproducibility-protocol methodology this project follows, three
levels are worth distinguishing, and this project claims different levels
for different things:

- **Exact (bit-for-bit) reproducibility:** claimed and tested for any
  single Monte Carlo run given the same `(master_seed, n_iterations)` and
  the same package versions (`tests/unit/test_monte_carlo.py::TestMonteCarloEngine::test_reproducibility_same_seed`,
  `tests/integration/test_scenario_runs.py::test_full_reproducibility_across_full_pipeline`).
- **Statistical reproducibility (within Monte Carlo error):** claimed
  across different seeds at adequate iteration count — i.e. summary
  statistics (mean, percentiles) should agree within the Monte Carlo
  standard error reported by `simulations/convergence.py`, not exactly.
- **NOT claimed:** bit-for-bit reproducibility across different NumPy/SciPy
  major versions (their internal random-number consumption patterns are not
  guaranteed stable across versions) or across different hardware/OS in
  principle (in practice, NumPy's `Generator`-based RNGs are algorithmic and
  should be stable, but this project has not tested cross-platform).

## Seed management strategy

Every Monte Carlo run uses **one master seed**, expanded into per-iteration,
per-purpose child seeds via `numpy.random.SeedSequence.spawn()`
(`simulations/monte_carlo.py:run_monte_carlo`):

```python
seed_seq = np.random.SeedSequence(master_seed)
child_seeds = seed_seq.spawn(n_iterations * 2)   # 2 independent streams per iteration
```

**Why two streams per iteration, not one:** iteration `i`'s outer-loop
parameter sampling and inner-loop stochastic simulation draw from
independent `Generator` instances (`child_seeds[2*i]` and
`child_seeds[2*i+1]`). This means, e.g., turning on `store_results=True`
(which changes nothing about which random numbers are consumed, only
whether results are kept in memory) cannot silently change the numerical
results — a common reproducibility bug class this design rules out by
construction.

**Never use `np.random.seed()` or the global `np.random` module anywhere in
this codebase.** Every stochastic function takes an explicit
`rng: np.random.Generator` argument. This is enforced by convention (code
review) rather than a runtime check — a static-analysis rule to enforce it
automatically is a reasonable extension (see `docs/roadmap.md`).

## Environment specification

- `pyproject.toml` — the authoritative dependency specification (version
  floors, not exact pins, appropriate for a library/research-tool rather
  than a deployed application).
- `requirements.txt` — exact versions this project was developed and tested
  against, for users who want a known-working environment rather than
  "whatever satisfies the floors today."
- `Dockerfile` — a fully specified container environment for maximum
  reproducibility, including the base Python version and OS.

## How to reproduce this project's own reported numbers

```bash
pip install -e ".[dev]"
pytest tests/ -q                                    # 108 tests should pass
pytest tests/ --cov=outbreak_simulator               # coverage report (~88%)
python scripts/generate_evidence_tables.py           # regenerate docs/evidence_tables.md
python examples/01_run_choir_outbreak.py             # reproduces the choir/measles example numbers in this doc set
```

All example scripts and the validation test suite use fixed, stated seeds
(`20260720` = this project's build date, or explicit small integers) —
grep for `seed=` in `examples/` and `tests/validation/` to find every one.

## Non-determinism sources this project is aware of and does NOT control for

- **Floating-point summation order** in NumPy's vectorized operations can
  in principle vary across BLAS backends/hardware in ways that produce
  last-bit differences — not expected to matter for this project's
  precision requirements, but noted for completeness per the
  reproducibility-protocol methodology's "document non-determinism sources"
  step.
- **matplotlib figure rendering** (font hinting, anti-aliasing) is not
  bit-reproducible across systems and is not claimed to be — only the
  underlying data feeding a plot is covered by the reproducibility
  guarantees above.
