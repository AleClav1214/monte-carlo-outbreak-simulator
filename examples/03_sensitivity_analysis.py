"""
Example 3: One-way (tornado) and global (PRCC) sensitivity analysis for the
university dormitory scenario.

Run with:  python examples/03_sensitivity_analysis.py
"""

from outbreak_simulator.data import get_pathogen, get_scenario
from outbreak_simulator.sensitivity import one_way_sensitivity, partial_rank_correlation
from outbreak_simulator.simulations import run_scenario
from outbreak_simulator.visualization import prcc_bar_chart, tornado_chart

if __name__ == "__main__":
    scenario = get_scenario("university_dormitory")
    pathogen = get_pathogen(scenario.pathogen_id)

    print("=== One-way sensitivity (tornado) ===")
    tornado_results = one_way_sensitivity(
        population_size=scenario.population.population_size,
        initial_cases=scenario.population.initial_cases,
        baseline_r0=pathogen.parameters["r0"].point_estimate,
        baseline_k=pathogen.parameters["k_dispersion"].point_estimate,
        contact_multiplier=scenario.population.contact_rate_multiplier,
        parameters_to_vary={"r0": pathogen.parameters["r0"], "k_dispersion": pathogen.parameters["k_dispersion"]},
        n_reps_per_point=500,
    )
    for r in tornado_results:
        print(f"  {r.parameter_name}: [{r.output_at_low:.3f}, {r.output_at_high:.3f}]  (swing={r.swing:.3f})")
    tornado_chart(tornado_results, title="University dormitory: one-way sensitivity").savefig("tornado.png", dpi=120)

    print("\n=== Global sensitivity (PRCC) ===")
    run_result = run_scenario("university_dormitory", n_iterations=5000, seed=1)
    independent_params = {k: v for k, v in run_result.mc_result.sampled_parameters.items() if k != "r_effective"}
    prcc_results = partial_rank_correlation(independent_params, run_result.mc_result.raw_attack_rates)
    for p in prcc_results:
        significance = "***" if p.p_value < 0.001 else ("*" if p.p_value < 0.05 else "n.s.")
        print(f"  {p.parameter_name}: PRCC={p.prcc:+.3f}  p={p.p_value:.4f} {significance}")
    prcc_bar_chart(prcc_results, title="University dormitory: global sensitivity (PRCC)").savefig("prcc.png", dpi=120)

    print("\nSaved tornado.png and prcc.png")
