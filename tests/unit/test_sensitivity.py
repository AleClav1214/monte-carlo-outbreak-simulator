"""Unit tests: sensitivity/ (one-way tornado analysis and global PRCC)."""

from __future__ import annotations

import numpy as np
import pytest

from outbreak_simulator.sensitivity.global_sensitivity import (
    leave_one_out_variance_contribution,
    partial_rank_correlation,
    scenario_robustness_analysis,
)
from outbreak_simulator.sensitivity.one_way import one_way_sensitivity, one_way_sensitivity_generic


class TestOneWaySensitivity:
    def test_returns_result_per_parameter(self, sars_cov_2):
        results = one_way_sensitivity(
            population_size=200, initial_cases=1,
            baseline_r0=sars_cov_2.parameters["r0"].point_estimate,
            baseline_k=sars_cov_2.parameters["k_dispersion"].point_estimate,
            contact_multiplier=1.5,
            parameters_to_vary={
                "r0": sars_cov_2.parameters["r0"], "k_dispersion": sars_cov_2.parameters["k_dispersion"]
            },
            n_reps_per_point=50,
        )
        assert len(results) == 2
        expected_names = {"Basic reproduction number (R0)", "Negative-binomial dispersion parameter (k)"}
        assert {r.parameter_name for r in results} == expected_names

    def test_sorted_by_swing_descending(self, sars_cov_2):
        results = one_way_sensitivity(
            population_size=200, initial_cases=1,
            baseline_r0=sars_cov_2.parameters["r0"].point_estimate,
            baseline_k=sars_cov_2.parameters["k_dispersion"].point_estimate,
            contact_multiplier=1.5,
            parameters_to_vary={
                "r0": sars_cov_2.parameters["r0"], "k_dispersion": sars_cov_2.parameters["k_dispersion"]
            },
            n_reps_per_point=50,
        )
        swings = [r.swing for r in results]
        assert swings == sorted(swings, reverse=True)

    def test_rejects_unsupported_parameter_name(self, sars_cov_2):
        with pytest.raises(ValueError):
            one_way_sensitivity(
                population_size=200, initial_cases=1, baseline_r0=2.5, baseline_k=0.5, contact_multiplier=1.0,
                parameters_to_vary={"secondary_attack_rate": sars_cov_2.parameters["secondary_attack_rate"]},
                n_reps_per_point=10,
            )

    def test_generic_one_way_works_with_custom_simulate_fn(self):
        def toy_sim(params: dict, seed: int) -> float:
            return params["x"] * 2  # deterministic, so sensitivity should exactly track x

        results = one_way_sensitivity_generic(
            baseline_params={"x": 5.0}, param_ranges={"x": (1.0, 10.0)},
            simulate_fn=toy_sim, n_reps_per_point=5,
        )
        assert results[0].output_at_low == pytest.approx(2.0)
        assert results[0].output_at_high == pytest.approx(20.0)


class TestPRCC:
    def test_requires_at_least_two_parameters(self):
        with pytest.raises(ValueError):
            partial_rank_correlation(
                {"x": np.random.default_rng(1).normal(size=100)}, np.random.default_rng(2).normal(size=100)
            )

    def test_detects_strong_known_driver(self):
        """A synthetic case where output is a deterministic (noisy) function of
        x1 only -- PRCC for x1 should be much stronger than for the irrelevant x2."""
        rng = np.random.default_rng(1)
        x1 = rng.uniform(size=2000)
        x2 = rng.uniform(size=2000)  # irrelevant
        output = x1 * 10 + rng.normal(scale=0.1, size=2000)
        results = partial_rank_correlation({"x1": x1, "x2": x2}, output)
        by_name = {r.parameter_name: r for r in results}
        assert abs(by_name["x1"].prcc) > 0.9
        assert abs(by_name["x2"].prcc) < 0.2

    def test_prcc_bounded(self):
        rng = np.random.default_rng(1)
        params = {"a": rng.normal(size=500), "b": rng.normal(size=500), "c": rng.normal(size=500)}
        output = params["a"] + rng.normal(scale=0.5, size=500)
        for r in partial_rank_correlation(params, output):
            assert -1.0 <= r.prcc <= 1.0


class TestLeaveOneOutAndScenarios:
    def test_leave_one_out_variance_contributions_sum_reasonably(self, sars_cov_2):
        from scipy import stats
        r0 = sars_cov_2.parameters["r0"]
        k = sars_cov_2.parameters["k_dispersion"]
        dists = {
            "r0": stats.gamma(
                a=(r0.point_estimate / ((r0.high - r0.low) / 4)) ** 2,
                scale=((r0.high - r0.low) / 4) ** 2 / r0.point_estimate,
            ),
            "k_dispersion": stats.lognorm(s=0.8, scale=k.point_estimate),
        }
        results = leave_one_out_variance_contribution(
            population_size=100, initial_cases=1,
            baseline_point_estimates={"r0": r0.point_estimate, "k_dispersion": k.point_estimate},
            parameter_distributions=dists, contact_multiplier=1.5, n_iterations=300,
        )
        assert len(results) == 2
        for r in results:
            # allow small negative from MC noise, but not wildly so
            assert -0.5 <= r.variance_contribution_fraction <= 1.0

    def test_scenario_robustness_orders_best_worse_case(self):
        results = scenario_robustness_analysis(
            population_size=100, initial_cases=1, contact_multiplier=1.0,
            scenarios={
                "best_case": {"r0": 1.5, "k_dispersion": 1.0},
                "worst_case": {"r0": 4.0, "k_dispersion": 1.0},
            },
            n_reps=300,
        )
        by_name = {r.scenario_name: r for r in results}
        assert by_name["worst_case"].mean_attack_rate > by_name["best_case"].mean_attack_rate
