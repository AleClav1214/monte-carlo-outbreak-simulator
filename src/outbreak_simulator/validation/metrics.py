"""
Goodness-of-fit / calibration metrics.

Standalone metric functions, kept separate from the calibration workflow
(calibration.py) that applies them to specific scenarios -- these are pure
functions operating on arrays, independently testable and reusable.

A note on what these metrics can and cannot tell you: with n=1 (or a
handful) of real-world observed outcomes per scenario, no metric here
constitutes "validation" in the sense of demonstrating the model
generalizes -- see docs/validation_plan.md for the internal/external/
functional validation-tier framework this project uses to avoid overclaiming.
What these metrics DO support: checking whether a specific, real, literature-
reported outcome is a *plausible* draw from the model's predictive
distribution, which is a necessary (not sufficient) condition for the model
being reasonable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CoverageResult:
    """Does the model's predictive interval contain the observed value?"""

    observed_value: float
    predictive_interval: tuple[float, float]
    interval_level: float  # e.g. 0.95
    covered: bool
    observed_percentile: float  # what percentile of the simulated distribution the observed value falls at


def predictive_coverage(simulated: np.ndarray, observed_value: float, interval_level: float = 0.95) -> CoverageResult:
    """Check whether `observed_value` falls within the simulated distribution's
    central interval_level interval, and report its percentile rank within
    the simulated distribution (useful even when the interval is very wide
    due to bimodality -- see docs/validation_plan.md discussion of the choir
    scenario's bimodal output)."""
    alpha = (1 - interval_level) / 2
    lo, hi = np.percentile(simulated, [100 * alpha, 100 * (1 - alpha)])
    covered = lo <= observed_value <= hi
    percentile = float((simulated <= observed_value).mean() * 100)
    return CoverageResult(
        observed_value=observed_value, predictive_interval=(float(lo), float(hi)),
        interval_level=interval_level, covered=bool(covered), observed_percentile=percentile,
    )


def root_mean_squared_error(simulated_point_estimates: np.ndarray, observed_values: np.ndarray) -> float:
    """RMSE between paired simulated point-estimates and observed values --
    meaningful when comparing multiple scenarios' mean/median predictions
    against multiple corresponding observed outcomes (not meaningful for a
    single scenario's single observed point, which is what coverage/bias
    below are for)."""
    sim = np.asarray(simulated_point_estimates, dtype=float)
    obs = np.asarray(observed_values, dtype=float)
    if sim.shape != obs.shape:
        raise ValueError("simulated and observed arrays must have matching shape")
    return float(np.sqrt(np.mean((sim - obs) ** 2)))


def relative_bias(simulated_point_estimate: float, observed_value: float) -> float:
    """(simulated - observed) / observed -- signed, so positive means the
    model over-predicts relative to this specific observation."""
    if observed_value == 0:
        return float("nan")
    return (simulated_point_estimate - observed_value) / observed_value


def posterior_predictive_check(simulated: np.ndarray, observed_value: float) -> float:
    """Two-sided posterior-predictive-check-style p-value: probability, under
    the model, of an outcome at least as extreme (in either direction) as
    the one observed. Very small values (e.g. <0.05) suggest the observed
    outcome is a poor fit to the model's predictive distribution; NOTE this
    is a diagnostic heuristic, not a formal hypothesis test with a
    pre-registered null -- do not report it as a p-value from a designed
    experiment."""
    sim = np.asarray(simulated, dtype=float)
    percentile = (sim <= observed_value).mean()
    return float(2 * min(percentile, 1 - percentile))
