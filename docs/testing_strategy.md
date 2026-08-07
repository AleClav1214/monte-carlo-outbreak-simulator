# Testing Strategy

## 1. Test organization

```
tests/
├── unit/          # one module's behavior in isolation, fast (<1s per test typically)
├── integration/   # multiple modules together: full scenario pipeline, intervention + MC + validation
└── validation/    # correctness against real literature benchmarks and known ground truth
```

**108 tests, 88% statement coverage** as of this writing (run
`pytest tests/ --cov=outbreak_simulator --cov-report=term-missing` to
reproduce). The uncovered 12% is concentrated almost entirely in
`transmission_inference/transphylo_interface.py` (33% covered) — see §4.

## 2. What "passing" means for each test tier

- **Unit tests** assert exact or statistically-precise properties: negative
  binomial offspring mean/variance matching closed-form theory to within
  Monte Carlo error, Beta/Gamma/Lognormal samples matching their target
  mean to within 1-5%, population conservation in the SEIR model being
  *exact* (not approximate), reproducibility being *bit-exact* under a
  fixed seed.
- **Integration tests** assert pipeline-level properties: every bundled
  scenario runs without error, interventions measurably reduce attack rate
  in the expected direction, PRCC runs successfully on real scenario output.
- **Validation tests** assert real-world consistency: observed benchmark
  attack rates are not statistical outliers relative to the model's
  predictive distribution (posterior-predictive-check p-value above a
  threshold), and the transmission-reconstruction method beats a
  random-guess baseline on synthetic ground truth at realistic scale.

## 3. Deliberate test design choices worth noting

- **`test_visualization.py` uses the `Agg` (headless) matplotlib backend**
  and asserts figures are produced without error — it does not do
  pixel-level image comparison. Visual correctness (e.g. "does the network
  plot correctly highlight the true index case as a superspreader") was
  checked manually during development (see commit history / development
  notes) rather than automated, since image-diff testing has a poor
  cost/benefit ratio for a project this size and tends to be brittle across
  matplotlib versions.
- **PRCC tests include a synthetic known-driver case**
  (`test_sensitivity.py::TestPRCC::test_detects_strong_known_driver`):
  construct output as a noisy deterministic function of one parameter and
  assert PRCC correctly identifies it (>0.9) while an irrelevant second
  parameter scores near zero (<0.2). This is what actually validates the
  PRCC *implementation* is correct, as opposed to merely "runs without
  crashing" on real (noisier, multi-cause) scenario data.
- **The `testing_isolation` naming collision.** An early version of
  `tests/unit/test_interventions.py` imported
  `interventions.behavioral.testing_isolation` directly, which pytest's
  default collection then mistook for a test function (name starts with
  `test`) and tried to call with no arguments — an immediate, confusing
  collection error. Fixed by aliasing the import
  (`testing_isolation as apply_testing_isolation`). Documented here because
  it is a genuinely easy mistake to reintroduce if a future intervention
  function is also named starting with `test_`/`Test`.
- **Reconstruction accuracy is asserted with a loose bound (>15%), not a
  tight one.** See `docs/validation_plan.md` §4 — the real, measured
  accuracy at realistic scale is ~25-38% across repeated trials, which is
  itself noisy (small-outbreak sample sizes). Asserting a tight bound would
  make the test suite flaky for the wrong reason (Monte Carlo noise in the
  *test's own* synthetic-outbreak generation) rather than catching genuine
  regressions.

## 4. Known coverage gaps and why

- **`transphylo_interface.py` (33% covered):** this module requires R +
  the TransPhylo package, which are not installed in this project's default
  environment (they are an optional extra — see `pyproject.toml`'s
  `transmission-inference-r` group) and are not available in this
  development sandbox. The covered 33% is the pure-Python logic that
  doesn't require R (config construction, the years-conversion helper, and
  the "gracefully raise a clear error when R isn't available" path, which
  IS tested — see `TransPhyloNotAvailableError` handling). The R-dependent
  code paths (`run_transphylo_rpy2`, `run_transphylo_subprocess`'s success
  path) are reviewed by inspection and documented with exact expected R
  function signatures, but not executed in CI. **This is stated here rather
  than hidden** — do not assume the R integration path is bug-free without
  independently testing it against a real TransPhylo installation.

## 5. Continuous Integration

`.github/workflows/ci.yml` runs the full test suite (excluding anything
requiring R) on Python 3.10, 3.11, and 3.12 on every push and pull request,
plus `ruff` linting. See that file directly for the exact matrix.

## 6. How to add a test for a new scenario or pathogen

1. Add the pathogen/scenario to the relevant YAML file — schema validation
   (`tests/unit/test_data_validation.py`) will automatically enforce the
   required fields.
2. `TestAllBundledScenariosRun` in `tests/integration/test_scenario_runs.py`
   automatically picks up any new scenario via `list_scenarios()` — no new
   test code needed for the basic "does it run" check.
3. If a real-world benchmark exists, add it to the scenario's
   `observed_outcomes` — `TestAllBundledScenariosRun::test_every_scenario_produces_a_calibration_report`
   will then enforce that at least one benchmark is present, and you should
   add a scenario-specific method to
   `tests/validation/test_calibration_against_literature.py` mirroring the
   existing ones.
