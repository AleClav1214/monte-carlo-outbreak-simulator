"""Unit tests: models/ (distributions, branching process, SEIR)."""

from __future__ import annotations

import numpy as np
import pytest

from outbreak_simulator.data import get_parameter
from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel
from outbreak_simulator.models.distributions import (
    beta_from_mean_ci,
    gamma_from_mean_ci,
    lognormal_from_mean_ci,
    negative_binomial_offspring,
    sample_parameter,
)
from outbreak_simulator.models.seir import SEIRConfig, StochasticSEIRModel


class TestDistributionFitting:
    def test_beta_from_mean_ci_matches_target_mean(self, rng):
        dist = beta_from_mean_ci(0.2, 0.1, 0.3)
        samples = dist.rvs(size=200_000, random_state=rng)
        assert abs(samples.mean() - 0.2) < 0.01

    def test_beta_from_mean_ci_bounded_unit_interval(self, rng):
        dist = beta_from_mean_ci(0.5, 0.3, 0.7)
        samples = dist.rvs(size=10_000, random_state=rng)
        assert samples.min() >= 0 and samples.max() <= 1

    def test_gamma_from_mean_ci_matches_target_mean(self, rng):
        dist = gamma_from_mean_ci(5.0, 3.0, 7.0)
        samples = dist.rvs(size=200_000, random_state=rng)
        assert abs(samples.mean() - 5.0) < 0.05

    def test_gamma_is_nonnegative(self, rng):
        dist = gamma_from_mean_ci(5.0, 3.0, 7.0)
        assert dist.rvs(size=10_000, random_state=rng).min() >= 0

    def test_lognormal_from_mean_ci_matches_target_mean(self, rng):
        dist = lognormal_from_mean_ci(5.8, 5.0, 6.7)
        samples = dist.rvs(size=200_000, random_state=rng)
        assert abs(samples.mean() - 5.8) < 0.05

    def test_sample_parameter_reproduces_evidence_table_mean(self, sars_cov_2, rng):
        r0_param = sars_cov_2.parameters["r0"]
        samples = sample_parameter(r0_param, rng, size=100_000)
        assert abs(samples.mean() - r0_param.point_estimate) < 0.02


class TestNegativeBinomialOffspring:
    def test_mean_matches_r(self, rng):
        offspring = negative_binomial_offspring(r_effective=2.5, k=0.5, rng=rng, size=200_000)
        assert abs(offspring.mean() - 2.5) < 0.05

    def test_variance_exceeds_poisson_variance_lowk(self, rng):
        """Overdispersion check: for k < infinity, Var > Mean (unlike Poisson, where Var == Mean)."""
        offspring = negative_binomial_offspring(r_effective=2.5, k=0.15, rng=rng, size=200_000)
        assert offspring.var() > offspring.mean() * 5  # dramatically overdispersed at k=0.15

    def test_rejects_nonpositive_k(self, rng):
        with pytest.raises(ValueError):
            negative_binomial_offspring(r_effective=2.5, k=0, rng=rng)

    def test_rejects_negative_r(self, rng):
        with pytest.raises(ValueError):
            negative_binomial_offspring(r_effective=-1.0, k=0.5, rng=rng)

    def test_zero_r_gives_all_zero_offspring(self, rng):
        offspring = negative_binomial_offspring(r_effective=0.0, k=0.5, rng=rng, size=1000)
        assert (offspring == 0).all()


class TestBranchingProcessModel:
    def test_rejects_invalid_population_size(self):
        cfg = BranchingProcessConfig(population_size=0, initial_cases=1, r_effective=2.0, k_dispersion=0.5)
        with pytest.raises(ValueError):
            BranchingProcessModel(cfg)

    def test_rejects_initial_cases_exceeding_population(self):
        cfg = BranchingProcessConfig(population_size=10, initial_cases=20, r_effective=2.0, k_dispersion=0.5)
        with pytest.raises(ValueError):
            BranchingProcessModel(cfg)

    def test_final_size_never_exceeds_population(self, rng):
        cfg = BranchingProcessConfig(population_size=50, initial_cases=1, r_effective=10.0, k_dispersion=1.0)
        model = BranchingProcessModel(cfg)
        for _ in range(200):
            result = model.run(rng)
            assert result.final_size <= 50
            assert 0 <= result.attack_rate <= 1.0

    def test_zero_r_effective_never_spreads_beyond_index_cases(self, rng):
        cfg = BranchingProcessConfig(population_size=100, initial_cases=3, r_effective=0.0, k_dispersion=0.5)
        model = BranchingProcessModel(cfg)
        result = model.run(rng)
        assert result.final_size == 3

    def test_very_high_r_in_small_population_saturates(self, rng):
        """A very high R in a small population should reliably infect nearly everyone (once it establishes)."""
        cfg = BranchingProcessConfig(population_size=30, initial_cases=1, r_effective=20.0, k_dispersion=5.0)
        model = BranchingProcessModel(cfg)
        rates = [model.run(rng).attack_rate for _ in range(100)]
        assert np.mean(rates) > 0.5  # high enough R/k that most realizations should take off and saturate

    def test_reproducibility_same_seed_same_result(self):
        cfg = BranchingProcessConfig(population_size=100, initial_cases=1, r_effective=3.0, k_dispersion=0.3)
        model = BranchingProcessModel(cfg)
        r1 = model.run(np.random.default_rng(42))
        r2 = model.run(np.random.default_rng(42))
        assert r1.final_size == r2.final_size
        assert r1.attack_rate == r2.attack_rate


class TestStochasticSEIRModel:
    def test_population_conserved(self, rng):
        cfg = SEIRConfig(population_size=200, initial_infectious=2, beta=0.4, sigma=0.2, gamma=0.15, max_time=100)
        model = StochasticSEIRModel(cfg)
        result = model.run(rng)
        final_state = result.metadata["final_state"]
        assert sum(final_state[c] for c in ("S", "E", "I", "R")) == 200

    def test_r0_zero_beta_causes_no_onward_spread(self, rng):
        cfg = SEIRConfig(population_size=200, initial_infectious=2, beta=0.0, sigma=0.2, gamma=0.15, max_time=100)
        model = StochasticSEIRModel(cfg)
        result = model.run(rng)
        # final_size includes the initial seed cases themselves (they are infected,
        # just not via within-model transmission); with beta=0 no ADDITIONAL S->E
        # transitions can occur, so final_size must equal exactly initial_infectious.
        assert result.final_size == 2
        assert result.metadata["final_state"]["S"] == 198  # no susceptible was ever converted

    def test_high_r0_produces_large_outbreak(self, rng):
        cfg = SEIRConfig(population_size=500, initial_infectious=3, beta=2.0, sigma=0.2, gamma=0.15, max_time=300)
        model = StochasticSEIRModel(cfg)
        result = model.run(rng)
        assert result.attack_rate > 0.5  # R0 = beta/gamma = 13.3, should produce a large epidemic

    def test_rejects_initial_infections_exceeding_population(self):
        cfg = SEIRConfig(population_size=10, initial_infectious=15)
        with pytest.raises(ValueError):
            StochasticSEIRModel(cfg)
