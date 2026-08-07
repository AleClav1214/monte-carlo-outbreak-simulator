# Contributing

## Development setup

```bash
git clone <this-repository>
cd monte-carlo-outbreak-simulator
pip install -e ".[dev]"
pytest tests/ -q   # should show 108 passed
```

## Before opening a pull request

1. **`pytest tests/`** must pass.
2. **`ruff check src/ tests/`** must pass (config in `pyproject.toml`).
3. If you changed `data/parameters/pathogens.yaml` or
   `data/scenarios/scenarios.yaml`, run
   `python scripts/generate_evidence_tables.py` and commit the regenerated
   `docs/evidence_tables.md`.
4. If you added a scenario or pathogen, confirm
   `tests/integration/test_scenario_runs.py::TestAllBundledScenariosRun`
   picks it up automatically (it should, via `list_scenarios()`) and add a
   real-world `observed_outcomes` entry if one exists in the literature.

## Contribution guidelines specific to this project

### Adding or changing a parameter estimate

Every `ParameterEstimate` requires `source`, `justification`, and
`evidence_quality` (see `data/schemas.py`). **Do not add a parameter
without a real, checkable source.** If you cannot find a dedicated
literature estimate, use `evidence_quality: low_confidence` and say so
explicitly in `justification` — this project's credibility rests on never
presenting a placeholder as if it were measured (see the several existing
`low_confidence` entries in `pathogens.yaml` for the expected tone and
level of detail).

### Adding a new intervention

Implement it as a factory function returning an `Intervention`
(`interventions/base.py`) in whichever of `pharmaceutical.py`,
`environmental.py`, or `behavioral.py` fits, with a real
`source` string. It will automatically compose with existing interventions
via `InterventionStack` — no changes needed there.

### Adding a new model

Implement the `OutbreakModel` interface (`models/base.py`: a `run(rng)`
method returning a `SimulationResult`, and a `population_size` property).
It will then work with `simulations/monte_carlo.py::run_monte_carlo`
without any changes to the Monte Carlo engine.

### Code style

- Every stochastic function takes an explicit `rng: np.random.Generator`
  parameter. Never call `np.random.*` module-level functions or construct
  `np.random.default_rng()` inside a function that's meant to be
  reproducible-by-caller — see `docs/reproducibility.md`.
- Prefer adding a test that would have failed before your fix, for any bug
  fix — several of this project's existing tests exist precisely because
  they caught a real bug during development (see `docs/testing_strategy.md`
  §3 for two examples: the `testing_isolation` pytest-collection naming
  collision, and a `numpy.bool_ is False` identity-comparison gotcha).

## Reporting issues

Use the issue templates in `.github/ISSUE_TEMPLATE/`. For a data/parameter
concern specifically ("this R0 estimate looks wrong" / "there's a newer
study"), please include the specific source you'd propose instead — this
project's evidence tables are only as good as what gets checked against
real literature.

## Code of conduct

Be respectful and constructive. This is a scientific tool; disagreements
about parameter values or modeling choices should be resolved by pointing
to evidence, not by assertion.
