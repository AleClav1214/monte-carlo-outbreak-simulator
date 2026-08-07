"""Monte Carlo simulation engine: parameter sampling, convergence diagnostics, and the high-level scenario runner."""

from outbreak_simulator.simulations.convergence import ConvergenceReport, assess_convergence, recommend_iteration_count
from outbreak_simulator.simulations.monte_carlo import (
    MonteCarloConfig,
    MonteCarloResult,
    OutputSummary,
    latin_hypercube_unit_samples,
    run_monte_carlo,
)
from outbreak_simulator.simulations.runner import ScenarioRunResult, print_summary, run_scenario

__all__ = [
    "ConvergenceReport", "assess_convergence", "recommend_iteration_count",
    "MonteCarloConfig", "MonteCarloResult", "OutputSummary", "latin_hypercube_unit_samples", "run_monte_carlo",
    "ScenarioRunResult", "print_summary", "run_scenario",
]
