"""
Monte Carlo uncertainty-propagation engine.

Two-level sampling design (see models/distributions.py for the full
rationale):
  - OUTER loop (this module): draws one value of each uncertain input
    parameter per iteration, from its evidence-table distribution --
    "what if the true R0 is 2.9 instead of 2.5?"
  - INNER loop (models/branching_process.py or models/seir.py): given that
    iteration's fixed parameter values, runs one stochastic realization of
    the transmission process -- "given R0=2.9, how many people actually get
    infected in this particular realization?"

This separates *epistemic* uncertainty (we don't know the true parameter
value) from *aleatory* stochasticity (even with known parameters, outcomes
are random) -- conflating the two is a common and consequential error in
outbreak modeling, and keeping them distinct is one of this project's
central design decisions.

Reproducibility: every iteration's randomness is derived from a single
master seed via numpy's SeedSequence.spawn(), which produces statistically
independent, non-overlapping child streams per iteration. This is the
modern numpy-recommended pattern (superseding global np.random.seed()) and
is what makes `scripts/reproduce.sh` (docs/reproducibility.md) able to
reproduce bit-identical results.

Output statistics: mean, median, variance, std, coefficient of variation,
skewness, kurtosis, and percentiles [5, 25, 50, 75, 95] plus the 95%
uncertainty interval -- exactly the set specified by this project's
Monte Carlo sampling methodology (docs/design_document.md).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from outbreak_simulator.models.base import OutbreakModel, SimulationResult


@dataclass
class MonteCarloConfig:
    n_iterations: int = 2000
    master_seed: int = 20260720  # date-stamped default; always pass an explicit seed for real analyses
    # if set, subsample SimulationResults kept in memory (summary stats use all iterations regardless)
    max_stored_results: int | None = None


@dataclass
class OutputSummary:
    """Required Monte Carlo output statistics for a single scalar output (e.g. attack_rate)."""

    mean: float
    median: float
    variance: float
    std: float
    coefficient_of_variation: float
    skewness: float
    kurtosis: float
    percentile_5: float
    percentile_25: float
    percentile_50: float
    percentile_75: float
    percentile_95: float
    uncertainty_interval_95: tuple[float, float]
    n: int

    @classmethod
    def from_samples(cls, samples: np.ndarray) -> OutputSummary:
        samples = np.asarray(samples, dtype=float)
        mean = float(np.mean(samples))
        std = float(np.std(samples, ddof=1)) if len(samples) > 1 else 0.0
        p = np.percentile(samples, [5, 25, 50, 75, 95])
        return cls(
            mean=mean,
            median=float(np.median(samples)),
            variance=float(np.var(samples, ddof=1)) if len(samples) > 1 else 0.0,
            std=std,
            coefficient_of_variation=(std / mean) if mean != 0 else float("nan"),
            skewness=float(stats.skew(samples)) if len(samples) > 2 else float("nan"),
            kurtosis=float(stats.kurtosis(samples)) if len(samples) > 3 else float("nan"),
            percentile_5=float(p[0]),
            percentile_25=float(p[1]),
            percentile_50=float(p[2]),
            percentile_75=float(p[3]),
            percentile_95=float(p[4]),
            uncertainty_interval_95=(float(p[0]), float(p[4])),
            n=len(samples),
        )

    def to_dict(self) -> dict:
        return {
            "mean": self.mean, "median": self.median, "variance": self.variance, "std": self.std,
            "coefficient_of_variation": self.coefficient_of_variation, "skewness": self.skewness,
            "kurtosis": self.kurtosis, "percentile_5": self.percentile_5, "percentile_25": self.percentile_25,
            "percentile_50": self.percentile_50, "percentile_75": self.percentile_75,
            "percentile_95": self.percentile_95, "uncertainty_interval_95": self.uncertainty_interval_95,
            "n": self.n,
        }


@dataclass
class MonteCarloResult:
    config: MonteCarloConfig
    attack_rate_summary: OutputSummary
    final_size_summary: OutputSummary
    peak_incidence_summary: OutputSummary
    duration_summary: OutputSummary
    extinction_probability: float
    raw_attack_rates: np.ndarray
    raw_final_sizes: np.ndarray
    stored_results: list[SimulationResult] = field(default_factory=list)
    # parameter name -> array of per-iteration draws
    sampled_parameters: dict[str, np.ndarray] = field(default_factory=dict)


ModelFactory = Callable[[np.random.Generator], tuple[OutbreakModel, dict]]
# A ModelFactory receives an rng (for drawing this iteration's OUTER-loop
# parameter values) and returns (constructed_model, sampled_params_dict).
# The returned model is then run with a *separate* child rng for its own
# INNER-loop stochasticity -- see run_monte_carlo() below for why these are
# deliberately different Generator instances.


def run_monte_carlo(
    model_factory: ModelFactory,
    config: MonteCarloConfig | None = None,
    store_results: bool = False,
) -> MonteCarloResult:
    """Run a full Monte Carlo batch.

    Two independent child RNG streams are spawned per iteration (one for
    outer-loop parameter sampling, one for inner-loop stochastic simulation)
    so that re-running with store_results=True/False, or re-running just the
    inner simulation with a fixed parameter draw, reproduces identical
    numbers -- the two concerns cannot accidentally consume each other's
    random draws.
    """
    cfg = config or MonteCarloConfig()
    seed_seq = np.random.SeedSequence(cfg.master_seed)
    child_seeds = seed_seq.spawn(cfg.n_iterations * 2)

    attack_rates = np.empty(cfg.n_iterations)
    final_sizes = np.empty(cfg.n_iterations)
    peak_incidences = np.empty(cfg.n_iterations)
    durations = np.empty(cfg.n_iterations)
    extinctions = np.empty(cfg.n_iterations, dtype=bool)
    stored: list[SimulationResult] = []
    sampled_params: dict[str, list[float]] = {}

    for i in range(cfg.n_iterations):
        param_rng = np.random.default_rng(child_seeds[2 * i])
        run_rng = np.random.default_rng(child_seeds[2 * i + 1])

        model, params_this_iter = model_factory(param_rng)
        for k, v in params_this_iter.items():
            sampled_params.setdefault(k, []).append(v)

        result = model.run(run_rng)
        attack_rates[i] = result.attack_rate
        final_sizes[i] = result.final_size
        peak_incidences[i] = result.peak_incidence
        durations[i] = result.duration
        extinctions[i] = result.extinct

        if store_results and (cfg.max_stored_results is None or i < cfg.max_stored_results):
            stored.append(result)

    return MonteCarloResult(
        config=cfg,
        attack_rate_summary=OutputSummary.from_samples(attack_rates),
        final_size_summary=OutputSummary.from_samples(final_sizes),
        peak_incidence_summary=OutputSummary.from_samples(peak_incidences),
        duration_summary=OutputSummary.from_samples(durations),
        extinction_probability=float(np.mean(extinctions)),
        raw_attack_rates=attack_rates,
        raw_final_sizes=final_sizes,
        stored_results=stored,
        sampled_parameters={k: np.array(v) for k, v in sampled_params.items()},
    )


def latin_hypercube_unit_samples(n_samples: int, n_dims: int, rng: np.random.Generator) -> np.ndarray:
    """Latin Hypercube samples on the unit hypercube [0,1)^n_dims.

    Used by sensitivity/ for more efficient parameter-space coverage than
    simple random sampling at a given sample size -- per the Monte Carlo
    sampling methodology's guidance to prefer LHS "for efficiency over
    simple random sampling." run_monte_carlo() itself uses simple random
    sampling (via each parameter's own distribution) for the headline
    uncertainty propagation, because independent random draws per parameter
    are simpler to reason about for the primary uncertainty intervals;
    LHS is used specifically where sample-efficient *coverage* of the joint
    parameter space matters more than independence, i.e. sensitivity
    analysis (sensitivity/multi_way.py) and convergence diagnostics.
    """
    result = np.empty((n_samples, n_dims))
    for d in range(n_dims):
        cut_points = (np.arange(n_samples) + rng.uniform(size=n_samples)) / n_samples
        result[:, d] = rng.permutation(cut_points)
    return result
