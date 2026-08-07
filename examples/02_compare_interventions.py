"""
Example 2: Compare multiple intervention strategies against a baseline for
the military barracks scenario, and visualize the comparison.

Run with:  python examples/02_compare_interventions.py
"""

from outbreak_simulator.interventions import (
    InterventionStack,
    masking,
    no_intervention,
    occupancy_reduction,
    testing_isolation,
    ventilation,
)
from outbreak_simulator.simulations import run_scenario
from outbreak_simulator.visualization import intervention_comparison_bar

if __name__ == "__main__":
    scenarios = {
        "No intervention": no_intervention(),
        "Masking only": InterventionStack("masking", [masking(coverage=0.8, evidence_basis="observational")]),
        "Ventilation only": InterventionStack("ventilation", [ventilation(baseline_ach=4.0, improved_ach=8.0)]),
        "Testing + isolation": InterventionStack("testing", [
            testing_isolation(test_sensitivity=0.85, testing_frequency_per_infectious_period=2)
        ]),
        "Full package": InterventionStack("full", [
            masking(coverage=0.8, evidence_basis="observational"),
            ventilation(baseline_ach=4.0, improved_ach=8.0),
            testing_isolation(test_sensitivity=0.85, testing_frequency_per_infectious_period=2),
            occupancy_reduction(occupancy_fraction=0.7),
        ]),
    }

    results = {}
    for name, stack in scenarios.items():
        run = run_scenario("military_barracks", intervention_stack=stack, n_iterations=5000, seed=1)
        results[name] = run.mc_result.attack_rate_summary.mean
        print(f"{name:25s}  mean attack rate = {run.mc_result.attack_rate_summary.mean:.1%}   "
              f"extinction prob = {run.mc_result.extinction_probability:.1%}")

    fig = intervention_comparison_bar(results, title="Military barracks: intervention comparison")
    fig.savefig("intervention_comparison.png", dpi=120)
    print("\nSaved intervention_comparison.png")
