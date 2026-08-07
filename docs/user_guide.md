# User Guide

## Installation

```bash
git clone <this-repository>
cd monte-carlo-outbreak-simulator
pip install -e .                              # core install
pip install -e ".[dev]"                       # + pytest, ruff (for development)
pip install -e ".[transmission-inference-r]"  # + rpy2 (only if you have R + TransPhylo installed)
```

Requires Python ≥3.10. See `docs/reproducibility.md` for exact dependency
versions this project was developed and tested against.

## Quickstart: run a bundled scenario

```python
from outbreak_simulator.simulations import run_scenario, print_summary

result = run_scenario("choir_rehearsal", n_iterations=5000, seed=20260720)
print(print_summary(result))
```

Available scenarios (see `docs/evidence_tables.md` for the full evidence
behind each):

| scenario_id | Setting |
|---|---|
| `choir_rehearsal` | Indoor choir rehearsal, SARS-CoV-2 |
| `military_barracks` | Basic training congregate housing, SARS-CoV-2 |
| `university_dormitory` | Residence hall, SARS-CoV-2 |
| `mpox_gathering` | Sexual-network-connected gathering, mpox |
| `influenza_outbreak` | Closed institutional setting, seasonal influenza |
| `norovirus_outbreak` | School/summer camp, norovirus |
| `measles_school_outbreak` | Under-vaccinated school/community, measles |
| `varicella_school_outbreak` | School/childcare, varicella |

## Checking what a scenario assumes before trusting its output

**Always read a scenario's `assumptions` and `limitations` before using its
output.** Every scenario has both:

```python
from outbreak_simulator.data import get_scenario

s = get_scenario("choir_rehearsal")
print(s.assumptions)
print(s.limitations)
```

For the measles and varicella school scenarios specifically, you **must**
set `susceptible_fraction` explicitly — population vaccination coverage is
a policy-relevant input this project deliberately does not hardcode:

```python
result = run_scenario("measles_school_outbreak", susceptible_fraction=0.08, n_iterations=5000)
```

## Modeling interventions

```python
from outbreak_simulator.interventions import InterventionStack, masking, ventilation, sars_cov_2_vaccination
from outbreak_simulator.simulations import run_scenario

stack = InterventionStack("layered protection", [
    masking(coverage=0.8, evidence_basis="observational"),   # or evidence_basis="rct" for the more conservative estimate
    ventilation(baseline_ach=1.0, improved_ach=6.0),
    sars_cov_2_vaccination(coverage=0.7, variant_era="delta"),
])
result = run_scenario("choir_rehearsal", intervention_stack=stack, n_iterations=5000)
```

Available intervention factory functions: `pharmaceutical.measles_vaccination`,
`pharmaceutical.sars_cov_2_vaccination`, `pharmaceutical.generic_vaccination`,
`environmental.masking`, `environmental.ventilation`, `environmental.air_filtration`,
`behavioral.testing_isolation`, `behavioral.quarantine_contacts`,
`behavioral.occupancy_reduction`.

## Comparing scenarios / interventions

```python
from outbreak_simulator.interventions import compare_scenarios, no_intervention

comparison = compare_scenarios(baseline_r_effective=8.5, stacks={
    "baseline": no_intervention(),
    "masks only": InterventionStack("masks", [masking(coverage=0.8)]),
    "full package": stack,
})
```

## Sensitivity analysis

```python
from outbreak_simulator.data import get_pathogen, get_scenario
from outbreak_simulator.sensitivity import one_way_sensitivity, partial_rank_correlation
from outbreak_simulator.visualization import tornado_chart

scenario = get_scenario("military_barracks")
pathogen = get_pathogen(scenario.pathogen_id)
tornado_results = one_way_sensitivity(
    population_size=scenario.population.population_size,
    initial_cases=scenario.population.initial_cases,
    baseline_r0=pathogen.parameters["r0"].point_estimate,
    baseline_k=pathogen.parameters["k_dispersion"].point_estimate,
    contact_multiplier=scenario.population.contact_rate_multiplier,
    parameters_to_vary={"r0": pathogen.parameters["r0"], "k_dispersion": pathogen.parameters["k_dispersion"]},
)
tornado_chart(tornado_results).savefig("tornado.png")
```

For global (PRCC) sensitivity, use the parameter samples a Monte Carlo run
already drew (`result.mc_result.sampled_parameters`) — see
`examples/03_sensitivity_analysis.py` for the full pattern, and read
`sensitivity/global_sensitivity.py`'s `partial_rank_correlation` docstring
**carefully** before including a composite/derived parameter (like
`r_effective`) alongside its own constituent factors — doing so silently
breaks the result (documented in detail in that docstring).

## Validating a scenario against its real-world benchmark

```python
from outbreak_simulator.validation import calibrate_scenario, print_calibration_report

report = calibrate_scenario(result)
print(print_calibration_report(report))
```

Read `docs/validation_plan.md` before interpreting this output — it is a
consistency check against a small number of real benchmarks, not a
certification that the model is "validated" in a strong sense.

## Transmission chain reconstruction

```python
from outbreak_simulator.transmission_inference import wallinga_teunis_reconstruction
from outbreak_simulator.visualization import plot_transmission_network

# case_times: array of symptom-onset times for your actual outbreak's cases
network = wallinga_teunis_reconstruction(case_times, generation_interval_shape=4.0, generation_interval_scale=1.3)
plot_transmission_network(network).savefig("network.png")
```

If you have real genomic sequence data (a dated phylogeny), see
`transmission_inference/transphylo_interface.py` instead — it requires R +
the TransPhylo package (`pip install ".[transmission-inference-r]"` plus R
setup, documented in that module's docstring) and is **not exercised by
this project's own examples or tests**, since no real genomic dataset is
bundled here (see `docs/validation_plan.md` §4 and
`docs/testing_strategy.md` §4 for why).

## Visualizing results

All plotting functions return a `matplotlib.figure.Figure` — call
`.savefig(path)` or embed in a notebook. See `visualization/__init__.py`
for the full list; `examples/` shows each in context.

## Extending this project with a new pathogen or scenario

1. Add a pathogen entry to `data/parameters/pathogens.yaml` with the five
   required parameters (`r0`, `incubation_period`, `infectious_period`,
   `k_dispersion`, `secondary_attack_rate`) — schema validation will reject
   the file if any are missing (`data/schemas.py`).
2. Add a scenario entry to `data/scenarios/scenarios.yaml` referencing that
   `pathogen_id` — referential integrity checks will reject a scenario
   pointing at an unknown pathogen.
3. Run `python scripts/generate_evidence_tables.py` to regenerate
   `docs/evidence_tables.md`.
4. Run `pytest tests/` — `TestAllBundledScenariosRun` automatically covers
   the new scenario.
