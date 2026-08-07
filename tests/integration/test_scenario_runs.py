"""Integration tests: full pipeline (data -> model -> interventions -> Monte Carlo -> validation)."""

from __future__ import annotations

from outbreak_simulator.data import list_scenarios
from outbreak_simulator.interventions import InterventionStack, masking, no_intervention, ventilation
from outbreak_simulator.sensitivity import partial_rank_correlation
from outbreak_simulator.simulations import run_scenario
from outbreak_simulator.validation import calibrate_scenario


class TestAllBundledScenariosRun:
    def test_every_scenario_runs_end_to_end(self):
        """Every bundled scenario must run without error at small iteration
        count -- this is the project's broadest regression test: if any
        scenario's data or parameter combination is broken, this catches it."""
        for scenario_id in list_scenarios():
            result = run_scenario(scenario_id, n_iterations=100, seed=1)
            assert result.mc_result.attack_rate_summary.n == 100
            assert 0 <= result.mc_result.attack_rate_summary.mean <= 1.0

    def test_every_scenario_produces_a_calibration_report(self):
        for scenario_id in list_scenarios():
            result = run_scenario(scenario_id, n_iterations=200, seed=1)
            report = calibrate_scenario(result)
            assert report.scenario_id == scenario_id
            # every bundled scenario has at least one observed_outcomes entry (data integrity check)
            assert report.n_independent_benchmarks >= 1, f"{scenario_id} has no external benchmark"


class TestInterventionIntegration:
    def test_interventions_reduce_attack_rate(self):
        baseline = run_scenario("choir_rehearsal", n_iterations=2000, seed=1, intervention_stack=no_intervention())
        stack = InterventionStack("protective measures", [
            masking(coverage=0.9, evidence_basis="observational"),
            ventilation(baseline_ach=1.0, improved_ach=6.0),
        ])
        protected = run_scenario("choir_rehearsal", n_iterations=2000, seed=1, intervention_stack=stack)
        assert protected.mc_result.attack_rate_summary.mean < baseline.mc_result.attack_rate_summary.mean
        assert protected.mc_result.extinction_probability > baseline.mc_result.extinction_probability

    def test_full_reproducibility_across_full_pipeline(self):
        r1 = run_scenario("military_barracks", n_iterations=500, seed=123)
        r2 = run_scenario("military_barracks", n_iterations=500, seed=123)
        assert r1.mc_result.attack_rate_summary.mean == r2.mc_result.attack_rate_summary.mean


class TestSensitivityIntegration:
    def test_prcc_runs_on_real_scenario_output(self):
        result = run_scenario("measles_school_outbreak", n_iterations=1000, seed=1)
        independent_params = {k: v for k, v in result.mc_result.sampled_parameters.items() if k != "r_effective"}
        prcc_results = partial_rank_correlation(independent_params, result.mc_result.raw_attack_rates)
        assert len(prcc_results) == len(independent_params)
        for p in prcc_results:
            assert -1.0 <= p.prcc <= 1.0
