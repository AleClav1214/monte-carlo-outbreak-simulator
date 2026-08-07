"""
Monte Carlo convergence diagnostics.

Answers the question "did I run enough iterations?" via three complementary
checks, per the Monte Carlo sampling methodology:
  1. Running-mean stabilization: does the cumulative mean stop drifting?
  2. Monte Carlo standard error (MCSE): the SEM of the estimate itself
     (std / sqrt(n)) -- shrinks as 1/sqrt(n), so this quantifies *how much*
     residual uncertainty is sampling noise vs. genuine parameter/model
     uncertainty.
  3. Split-chain R-hat-style diagnostic: split the run into two halves,
     compare between-half and within-half variance (a simplified analogue
     of the Gelman-Rubin R-hat used for MCMC convergence, adapted here for
     independent Monte Carlo draws rather than a Markov chain -- for truly
     independent draws this should already be close to 1.0, so a large
     deviation flags either a bug (e.g. accidental seed correlation) or a
     genuinely heavy-tailed/bimodal output distribution where "converged"
     needs a larger n than it would for a well-behaved unimodal output).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ConvergenceReport:
    n: int
    running_mean_final: float
    # how much the running mean moved over the final 10% of iterations, relative to its value
    running_mean_relative_change_last_10pct: float
    monte_carlo_standard_error: float
    split_half_rhat: float
    converged: bool
    recommendation: str


def running_mean(samples: np.ndarray) -> np.ndarray:
    return np.cumsum(samples) / np.arange(1, len(samples) + 1)


def monte_carlo_standard_error(samples: np.ndarray) -> float:
    return float(np.std(samples, ddof=1) / np.sqrt(len(samples)))


def split_half_rhat(samples: np.ndarray) -> float:
    """Simplified split-chain R-hat: split into two halves, compare between-half
    vs within-half variance. ~1.0 indicates no detectable inhomogeneity between
    the first and second half of the run (a necessary, not sufficient,
    condition for convergence)."""
    n = len(samples)
    half = n // 2
    a, b = samples[:half], samples[half : 2 * half]
    chain_means = np.array([a.mean(), b.mean()])
    chain_vars = np.array([a.var(ddof=1), b.var(ddof=1)])
    within = chain_vars.mean()
    between = half * chain_means.var(ddof=1)
    if within <= 0:
        return 1.0
    var_hat = ((half - 1) / half) * within + between / half
    return float(np.sqrt(var_hat / within))


def assess_convergence(
    samples: np.ndarray,
    relative_change_tolerance: float = 0.02,
    rhat_tolerance: float = 1.05,
) -> ConvergenceReport:
    samples = np.asarray(samples, dtype=float)
    n = len(samples)
    rm = running_mean(samples)
    tail_start = max(int(n * 0.9), 1)
    tail_change = abs(rm[-1] - rm[tail_start]) / (abs(rm[-1]) + 1e-12)
    mcse = monte_carlo_standard_error(samples)
    rhat = split_half_rhat(samples)

    converged = bool(tail_change < relative_change_tolerance and rhat < rhat_tolerance)
    if converged:
        recommendation = f"Converged: running mean moved {tail_change:.2%} over the final 10% of iterations (n={n})."
    elif rhat >= rhat_tolerance:
        recommendation = (
            f"NOT reliably converged: split-half R-hat={rhat:.3f} >= {rhat_tolerance}. This often indicates "
            f"a heavy-tailed or bimodal output (common for low-k superspreading scenarios on small populations "
            f"-- see docs/validation_plan.md) rather than a bug; increase n_iterations substantially and re-check "
            f"rather than trusting point estimates from this run."
        )
    else:
        recommendation = (
            f"NOT converged: running mean still moved {tail_change:.2%} over the final 10% of iterations "
            f"(tolerance {relative_change_tolerance:.2%}). Increase n_iterations."
        )
    return ConvergenceReport(
        n=n,
        running_mean_final=float(rm[-1]),
        running_mean_relative_change_last_10pct=float(tail_change),
        monte_carlo_standard_error=mcse,
        split_half_rhat=rhat,
        converged=converged,
        recommendation=recommendation,
    )


def recommend_iteration_count(
    pilot_samples: np.ndarray, target_mcse_relative: float = 0.01
) -> int:
    """Given a pilot Monte Carlo run, recommend an n_iterations that would
    bring the Monte Carlo standard error down to target_mcse_relative * mean
    (since MCSE shrinks as 1/sqrt(n), this scales the pilot's MCSE
    quadratically)."""
    pilot_samples = np.asarray(pilot_samples, dtype=float)
    n_pilot = len(pilot_samples)
    mean = np.mean(pilot_samples)
    if mean == 0:
        return n_pilot
    current_mcse_rel = monte_carlo_standard_error(pilot_samples) / abs(mean)
    if current_mcse_rel <= target_mcse_relative:
        return n_pilot
    scale = (current_mcse_rel / target_mcse_relative) ** 2
    return int(np.ceil(n_pilot * scale))
