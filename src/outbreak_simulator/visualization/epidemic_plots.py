"""
Epidemic curve and parameter-distribution visualization.

All functions return a matplotlib Figure (never call plt.show()) so callers
can save, embed, or further customize the figure -- standard library
convention, and necessary for non-interactive use (CI, batch report
generation).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from outbreak_simulator.data.schemas import ParameterEstimate
from outbreak_simulator.simulations.monte_carlo import MonteCarloResult


def epidemic_curve_with_uncertainty(
    mc_result: MonteCarloResult,
    title: str = "Simulated epidemic curve",
    max_time_bins: int = 60,
) -> plt.Figure:
    """Plot median daily incidence with a shaded 50%/95% uncertainty ribbon,
    built from the daily_incidence arrays of every stored SimulationResult
    in mc_result.stored_results (requires the Monte Carlo run to have been
    called with store_results=True)."""
    if not mc_result.stored_results:
        raise ValueError(
            "epidemic_curve_with_uncertainty requires per-iteration results; "
            "re-run run_monte_carlo(..., store_results=True)"
        )

    max_len = min(max(len(r.daily_incidence) for r in mc_result.stored_results), max_time_bins)
    padded = np.zeros((len(mc_result.stored_results), max_len))
    for i, r in enumerate(mc_result.stored_results):
        n = min(len(r.daily_incidence), max_len)
        padded[i, :n] = r.daily_incidence[:n]

    p5, p25, p50, p75, p95 = np.percentile(padded, [5, 25, 50, 75, 95], axis=0)
    x = np.arange(max_len)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(x, p5, p95, alpha=0.2, color="#c0392b", label="95% uncertainty interval")
    ax.fill_between(x, p25, p75, alpha=0.35, color="#c0392b", label="50% uncertainty interval")
    ax.plot(x, p50, color="#c0392b", linewidth=2, label="Median")
    ax.set_xlabel("Time (generation index or day, depending on model)")
    ax.set_ylabel("New infections")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def parameter_distribution_plot(param: ParameterEstimate, samples: np.ndarray | None = None) -> plt.Figure:
    """Plot a parameter's evidence-table interval alongside (optionally) the
    actual Monte Carlo samples drawn from it, so users can visually confirm
    the sampling distribution matches the literature-reported interval."""
    fig, ax = plt.subplots(figsize=(7, 4))
    if samples is not None:
        ax.hist(samples, bins=50, density=True, alpha=0.6, color="#2980b9", label="Monte Carlo samples")
    ax.axvline(param.point_estimate, color="#2c3e50", linestyle="-", linewidth=2, label="Point estimate")
    if param.low is not None and param.high is not None:
        ci_pct = int((param.ci_level or 0.95) * 100)
        ax.axvspan(param.low, param.high, alpha=0.15, color="#2c3e50", label=f"{ci_pct}% interval")
    ax.set_xlabel(f"{param.display_name} ({param.unit})")
    ax.set_ylabel("Density")
    ax.set_title(f"{param.display_name}\nsource: {param.source[:70]}{'...' if len(param.source) > 70 else ''}")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def attack_rate_histogram(
    mc_result: MonteCarloResult, observed_value: float | None = None, title: str = "Attack rate distribution"
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(mc_result.raw_attack_rates, bins=40, color="#8e44ad", alpha=0.7)
    s = mc_result.attack_rate_summary
    ax.axvline(s.median, color="#2c3e50", linewidth=2, label=f"Median ({s.median:.1%})")
    ax.axvspan(s.percentile_5, s.percentile_95, alpha=0.1, color="#2c3e50", label="95% UI")
    if observed_value is not None:
        ax.axvline(
            observed_value, color="#c0392b", linewidth=2, linestyle="--", label=f"Observed ({observed_value:.1%})"
        )
    ax.set_xlabel("Attack rate")
    ax.set_ylabel("Simulated iterations")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig
