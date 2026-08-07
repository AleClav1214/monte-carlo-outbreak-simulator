"""Unit tests: simulations/ (Monte Carlo engine, convergence diagnostics)."""

from __future__ import annotations

import numpy as np
import pytest

from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel
from outbreak_simulator.simulations.convergence import (
    assess_convergence,
    monte_carlo_standard_error,
    recommend_iteration_count,
    running_mean,
    split_half_rhat,
)
from outbreak_simulator.simulations.monte_carlo import (
    MonteCarloConfig,
    OutputSummary,
    latin_hypercube_unit_samples,
    run_monte_carlo,
)


def _simple_factory(rng):
    cfg = BranchingProcessConfig(population_size=100, initial_cases=1, r_effective=1.5, k_dispersion=0.5)
    return BranchingProcessModel(cfg), {"r_effective": 1.5}


class TestOutputSummary:
    def test_summary_statistics_correct_for_known_samples(self):
        samples = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        s = OutputSummary.from_samples(samples)
        assert s.mean == 3.0
        assert s.median == 3.0
        assert s.n == 5

    def test_all_required_statistics_present(self):
        """Regression test: the Monte Carlo skill's methodology requires
        mean, median, variance, std, CV, skewness, kurtosis, and the
        5/25/50/75/95 percentiles -- verify every field is populated."""
        s = OutputSummary.from_samples(np.random.default_rng(1).normal(size=1000))
        for field in ["mean", "median", "variance", "std", "coefficient_of_variation",
                      "skewness", "kurtosis", "percentile_5", "percentile_25",
                      "percentile_50", "percentile_75", "percentile_95", "uncertainty_interval_95"]:
            assert getattr(s, field) is not None

    def test_to_dict_roundtrips_keys(self):
        s = OutputSummary.from_samples(np.array([1.0, 2.0, 3.0]))
        d = s.to_dict()
        assert "mean" in d and "uncertainty_interval_95" in d


class TestMonteCarloEngine:
    def test_reproducibility_same_seed(self):
        r1 = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=500, master_seed=99))
        r2 = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=500, master_seed=99))
        assert r1.attack_rate_summary.mean == r2.attack_rate_summary.mean
        np.testing.assert_array_equal(r1.raw_attack_rates, r2.raw_attack_rates)

    def test_different_seeds_give_different_results(self):
        r1 = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=500, master_seed=1))
        r2 = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=500, master_seed=2))
        assert not np.array_equal(r1.raw_attack_rates, r2.raw_attack_rates)

    def test_n_iterations_respected(self):
        result = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=317, master_seed=1))
        assert len(result.raw_attack_rates) == 317

    def test_store_results_flag(self):
        result_stored = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=50, master_seed=1), store_results=True)
        result_not_stored = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=50, master_seed=1), store_results=False)
        assert len(result_stored.stored_results) == 50
        assert len(result_not_stored.stored_results) == 0

    def test_sampled_parameters_captured(self):
        result = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=10, master_seed=1))
        assert "r_effective" in result.sampled_parameters
        assert len(result.sampled_parameters["r_effective"]) == 10

    def test_extinction_probability_in_valid_range(self):
        result = run_monte_carlo(_simple_factory, MonteCarloConfig(n_iterations=500, master_seed=1))
        assert 0.0 <= result.extinction_probability <= 1.0


class TestConvergenceDiagnostics:
    def test_running_mean_converges_to_true_mean(self):
        rng = np.random.default_rng(1)
        samples = rng.normal(loc=5.0, scale=1.0, size=100_000)
        rm = running_mean(samples)
        assert abs(rm[-1] - 5.0) < 0.05

    def test_mcse_shrinks_with_more_samples(self):
        rng = np.random.default_rng(1)
        small = rng.normal(size=100)
        large = rng.normal(size=10_000)
        assert monte_carlo_standard_error(large) < monte_carlo_standard_error(small)

    def test_split_half_rhat_near_one_for_iid_samples(self):
        rng = np.random.default_rng(1)
        samples = rng.normal(size=50_000)
        rhat = split_half_rhat(samples)
        assert 0.95 < rhat < 1.1  # should be close to 1 for genuinely IID draws

    def test_assess_convergence_flags_small_n_as_not_converged(self):
        rng = np.random.default_rng(1)
        samples = rng.exponential(size=20)  # too few samples, heavy-tailed distribution
        report = assess_convergence(samples, relative_change_tolerance=0.001)
        assert report.converged is False

    def test_assess_convergence_flags_large_n_well_behaved_as_converged(self):
        rng = np.random.default_rng(1)
        samples = rng.normal(loc=10, scale=0.1, size=50_000)
        report = assess_convergence(samples)
        assert report.converged is True

    def test_recommend_iteration_count_increases_for_tighter_target(self):
        rng = np.random.default_rng(1)
        pilot = rng.normal(loc=10, scale=5, size=200)
        loose = recommend_iteration_count(pilot, target_mcse_relative=0.05)
        tight = recommend_iteration_count(pilot, target_mcse_relative=0.01)
        assert tight > loose


class TestLatinHypercubeSampling:
    def test_correct_shape(self):
        rng = np.random.default_rng(1)
        samples = latin_hypercube_unit_samples(100, 3, rng)
        assert samples.shape == (100, 3)

    def test_values_in_unit_interval(self):
        rng = np.random.default_rng(1)
        samples = latin_hypercube_unit_samples(50, 2, rng)
        assert samples.min() >= 0 and samples.max() < 1

    def test_stratification_one_point_per_bin(self):
        """Core LHS property: with n samples, each of the n equal-width bins
        along each dimension contains exactly one sample."""
        rng = np.random.default_rng(1)
        n = 20
        samples = latin_hypercube_unit_samples(n, 1, rng)
        bins = np.floor(samples[:, 0] * n).astype(int)
        assert sorted(bins) == list(range(n))
