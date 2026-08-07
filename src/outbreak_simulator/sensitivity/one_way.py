"""
One-way sensitivity analysis.

Varies ONE parameter at a time across its literature-reported uncertainty
range (5th/95th percentile of its distribution by default), holding all
other parameters at their point estimates, and measures the effect on a
chosen output statistic (default: mean attack rate). This is the standard
"tornado chart" methodology: rank parameters by the width of their output
swing to identify which uncertain inputs matter most for conclusions.

Deliberately simple and fast (a small number of point-estimate simulations
per parameter, not a full Monte Carlo batch per point) -- multi_way.py
covers the case where joint/interaction effects and full uncertainty
propagation across all parameters simultaneously matter.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from outbreak_simulator.data.schemas import ParameterEstimate
from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel


@dataclass
class OneWayResult:
    parameter_name: str
    low_value: float
    baseline_value: float
    high_value: float
    output_at_low: float
    output_at_baseline: float
    output_at_high: float

    @property
    def swing(self) -> float:
        """Total output range induced by this parameter's uncertainty -- the tornado bar width."""
        return abs(self.output_at_high - self.output_at_low)


def _mean_attack_rate_at_point(
    population_size: int, initial_cases: int, r_effective: float, k_dispersion: float,
    n_reps: int, seed: int,
) -> float:
    """Point-estimate output: mean attack rate across n_reps stochastic
    realizations at FIXED parameter values (this is inner-loop stochasticity
    only, no outer-loop parameter resampling -- appropriate for one-way
    sensitivity, where we want the output as a function of the parameter
    value, not blurred by simultaneously resampling everything else)."""
    rng = np.random.default_rng(seed)
    rates = np.empty(n_reps)
    for i in range(n_reps):
        cfg = BranchingProcessConfig(
            population_size=population_size, initial_cases=initial_cases,
            r_effective=r_effective, k_dispersion=k_dispersion,
        )
        rates[i] = BranchingProcessModel(cfg).run(rng).attack_rate
    return float(rates.mean())


def one_way_sensitivity(
    population_size: int,
    initial_cases: int,
    baseline_r0: float,
    baseline_k: float,
    contact_multiplier: float,
    parameters_to_vary: dict[str, ParameterEstimate],
    n_reps_per_point: int = 300,
    seed: int = 20260720,
) -> list[OneWayResult]:
    """Run one-way sensitivity analysis on r0 and k_dispersion (and any other
    named ParameterEstimates in parameters_to_vary -- keys must be 'r0' or
    'k_dispersion' since those are the two parameters the branching process
    model consumes directly).

    Returns results sorted by swing, descending -- ready to hand directly to
    visualization/sensitivity_plots.py:tornado_chart().
    """
    results = []
    for name, param in parameters_to_vary.items():
        if name not in ("r0", "k_dispersion"):
            raise ValueError(f"one_way_sensitivity only supports varying 'r0' or 'k_dispersion', got '{name}'")
        low = param.low if param.low is not None else param.point_estimate * 0.7
        high = param.high if param.high is not None else param.point_estimate * 1.3

        def output_for(r0_val: float, k_val: float, s: int) -> float:
            r_eff = r0_val * contact_multiplier
            return _mean_attack_rate_at_point(population_size, initial_cases, r_eff, k_val, n_reps_per_point, s)

        if name == "r0":
            out_low = output_for(low, baseline_k, seed)
            out_base = output_for(param.point_estimate, baseline_k, seed + 1)
            out_high = output_for(high, baseline_k, seed + 2)
        else:  # k_dispersion
            out_low = output_for(baseline_r0, low, seed)
            out_base = output_for(baseline_r0, baseline_k, seed + 1)
            out_high = output_for(baseline_r0, high, seed + 2)

        results.append(OneWayResult(
            parameter_name=param.display_name, low_value=low,
            baseline_value=param.point_estimate, high_value=high,
            output_at_low=out_low, output_at_baseline=out_base, output_at_high=out_high,
        ))

    return sorted(results, key=lambda r: r.swing, reverse=True)


def one_way_sensitivity_generic(
    baseline_params: dict[str, float],
    param_ranges: dict[str, tuple[float, float]],
    simulate_fn: Callable[[dict[str, float], int], float],
    n_reps_per_point: int = 300,
    seed: int = 20260720,
) -> list[OneWayResult]:
    """Generalized one-way sensitivity analysis for any simulate_fn(params_dict, seed) -> scalar output,
    not just the branching process model -- used e.g. for SEIR-model or
    intervention-parameter sensitivity where the specific dataclass config
    differs from BranchingProcessConfig."""
    results = []
    for name, (low, high) in param_ranges.items():
        baseline_val = baseline_params[name]

        def run_at(value: float, s: int, param_name: str = name) -> float:
            params = dict(baseline_params)
            params[param_name] = value
            outs = [simulate_fn(params, s + r) for r in range(n_reps_per_point)]
            return float(np.mean(outs))

        out_low = run_at(low, seed)
        out_base = run_at(baseline_val, seed + n_reps_per_point)
        out_high = run_at(high, seed + 2 * n_reps_per_point)
        results.append(OneWayResult(
            parameter_name=name, low_value=low, baseline_value=baseline_val, high_value=high,
            output_at_low=out_low, output_at_baseline=out_base, output_at_high=out_high,
        ))
    return sorted(results, key=lambda r: r.swing, reverse=True)
