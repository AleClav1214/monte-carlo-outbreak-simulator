"""Unit tests: visualization/ (smoke tests -- correctness of *content* is checked
visually during development; these tests guard against regressions that would
make plotting crash or silently produce empty figures)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless backend for CI

import numpy as np
import pytest

from outbreak_simulator.data import get_parameter
from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel
from outbreak_simulator.models.distributions import sample_parameter
from outbreak_simulator.sensitivity.global_sensitivity import PRCCResult
from outbreak_simulator.sensitivity.one_way import OneWayResult
from outbreak_simulator.simulations.monte_carlo import MonteCarloConfig, run_monte_carlo
from outbreak_simulator.transmission_inference import simulate_ground_truth_tree, wallinga_teunis_reconstruction
from outbreak_simulator.visualization import (
    attack_rate_histogram,
    epidemic_curve_with_uncertainty,
    intervention_comparison_bar,
    parameter_distribution_plot,
    plot_transmission_network,
    prcc_bar_chart,
    tornado_chart,
)


def _small_mc_result(store=True):
    def factory(rng):
        cfg = BranchingProcessConfig(population_size=60, initial_cases=1, r_effective=3.0, k_dispersion=0.3)
        return BranchingProcessModel(cfg), {"r_effective": 3.0}
    return run_monte_carlo(factory, MonteCarloConfig(n_iterations=100, master_seed=1), store_results=store)


class TestEpidemicPlots:
    def test_epidemic_curve_produces_figure(self):
        fig = epidemic_curve_with_uncertainty(_small_mc_result(store=True))
        assert fig is not None
        assert len(fig.axes) == 1

    def test_epidemic_curve_requires_stored_results(self):
        with pytest.raises(ValueError):
            epidemic_curve_with_uncertainty(_small_mc_result(store=False))

    def test_attack_rate_histogram_produces_figure(self):
        fig = attack_rate_histogram(_small_mc_result(), observed_value=0.5)
        assert fig is not None

    def test_parameter_distribution_plot_with_samples(self):
        param = get_parameter("sars_cov_2", "r0")
        samples = sample_parameter(param, np.random.default_rng(1), size=1000)
        fig = parameter_distribution_plot(param, samples)
        assert fig is not None

    def test_parameter_distribution_plot_without_samples(self):
        param = get_parameter("measles", "r0")
        fig = parameter_distribution_plot(param)
        assert fig is not None


class TestComparisonAndSensitivityPlots:
    def test_tornado_chart_produces_figure(self):
        results = [
            OneWayResult("R0", 2.0, 2.5, 3.0, 0.3, 0.4, 0.5),
            OneWayResult("k", 0.1, 0.15, 0.6, 0.2, 0.4, 0.9),
        ]
        fig = tornado_chart(results)
        assert fig is not None
        assert len(fig.axes[0].patches) == 2  # one bar per parameter

    def test_intervention_comparison_bar_produces_figure(self):
        fig = intervention_comparison_bar({"baseline": 0.5, "masks": 0.3, "masks+vent": 0.1})
        assert fig is not None

    def test_prcc_bar_chart_produces_figure(self):
        results = [PRCCResult("r0", 0.1, 0.3), PRCCResult("k", 0.6, 0.001)]
        fig = prcc_bar_chart(results)
        assert fig is not None


class TestNetworkPlots:
    def test_plot_transmission_network_nontrivial_outbreak(self):
        rng = np.random.default_rng(10)
        case_times, _ = simulate_ground_truth_tree(2.2, 0.5, 4.0, 1.3, 60, rng)
        network = wallinga_teunis_reconstruction(case_times, 4.0, 1.3)
        fig = plot_transmission_network(network)
        assert fig is not None

    def test_plot_transmission_network_handles_single_case(self):
        rng = np.random.default_rng(999)
        case_times, _ = simulate_ground_truth_tree(0.0, 0.5, 4.0, 1.3, 60, rng)  # r_eff=0 -> only index case
        network = wallinga_teunis_reconstruction(case_times, 4.0, 1.3)
        fig = plot_transmission_network(network)  # must not crash on a trivial 1-node network
        assert fig is not None
