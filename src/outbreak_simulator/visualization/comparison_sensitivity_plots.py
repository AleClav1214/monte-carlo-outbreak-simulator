"""Intervention-comparison and sensitivity-analysis (tornado chart) plots."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from outbreak_simulator.sensitivity.one_way import OneWayResult


def tornado_chart(results: list[OneWayResult], title: str = "One-way sensitivity analysis") -> plt.Figure:
    """Standard tornado chart: horizontal bars sorted by influence (swing),
    showing the output range induced by each parameter's uncertainty."""
    results = sorted(results, key=lambda r: r.swing)  # ascending, so largest ends up at top when plotted
    names = [r.parameter_name for r in results]
    lows = np.array([min(r.output_at_low, r.output_at_high) for r in results])
    highs = np.array([max(r.output_at_low, r.output_at_high) for r in results])
    baseline = results[0].output_at_baseline if results else 0.0

    fig, ax = plt.subplots(figsize=(8, max(3, 0.6 * len(results) + 1)))
    y = np.arange(len(results))
    ax.barh(y, highs - lows, left=lows, color="#2980b9", alpha=0.75, height=0.6)
    ax.axvline(baseline, color="#2c3e50", linestyle="--", linewidth=1.5, label="Baseline (all params at point estimate)")
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlabel("Output (mean attack rate)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def intervention_comparison_bar(scenario_results: dict[str, float], title: str = "Intervention scenario comparison", ylabel: str = "Mean attack rate") -> plt.Figure:
    """Bar chart comparing a metric (e.g. mean attack rate, or effective R)
    across named intervention scenarios -- the standard requirement-#6
    'scenario comparison' visualization."""
    names = list(scenario_results.keys())
    values = list(scenario_results.values())
    colors = plt.cm.RdYlGn_r(np.linspace(0.15, 0.85, len(names)))  # red=high risk, green=low risk

    fig, ax = plt.subplots(figsize=(max(6, 1.3 * len(names)), 5))
    bars = ax.bar(names, values, color=colors)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v:.2f}" if v > 1 else f"{v:.1%}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    ha="center", va="bottom", fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def prcc_bar_chart(prcc_results, title: str = "Global sensitivity (PRCC)") -> plt.Figure:
    """Bar chart of Partial Rank Correlation Coefficients, signed and sorted by magnitude."""
    prcc_results = sorted(prcc_results, key=lambda r: abs(r.prcc))
    names = [r.parameter_name for r in prcc_results]
    values = [r.prcc for r in prcc_results]
    colors = ["#c0392b" if v < 0 else "#27ae60" for v in values]

    fig, ax = plt.subplots(figsize=(7, max(3, 0.5 * len(names) + 1)))
    ax.barh(names, values, color=colors, alpha=0.8)
    ax.axvline(0, color="#2c3e50", linewidth=1)
    ax.set_xlabel("Partial Rank Correlation Coefficient")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig
