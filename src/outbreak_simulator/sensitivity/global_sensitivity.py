"""
Global (multi-way) sensitivity analysis.

Implements three complementary techniques, all operating on the joint
parameter samples already drawn by a Monte Carlo run (simulations/monte_carlo.py)
-- no additional simulation is required beyond what run_monte_carlo() already
produced, which is the efficiency argument for LHS/PRCC-based global
sensitivity analysis over a full factorial design:

  1. PRCC (Partial Rank Correlation Coefficient) -- the standard global
     sensitivity measure for nonlinear, monotonic input-output relationships
     in Monte Carlo epidemiological models (Marino, Hogue, Ray & Kirschner
     2008, J Theor Biol, "A methodology for performing global uncertainty
     and sensitivity analysis in systems biology" -- the de facto standard
     reference for this technique in exactly this application area). PRCC
     partials out the linear effect of every OTHER sampled parameter before
     computing each parameter's rank correlation with the output, isolating
     that parameter's own contribution even when parameters are examined
     jointly (unlike one-way analysis, which holds everything else fixed at
     a single point).

  2. Leave-one-out variance contribution -- for each parameter, re-run the
     Monte Carlo output-variance calculation with that parameter clamped to
     its point estimate (removing its contribution to output variance) and
     report the resulting *reduction* in output variance. This directly
     answers "how much of our total uncertainty comes from not knowing this
     parameter precisely?" -- distinct from PRCC, which answers "how
     strongly does this parameter drive the output level," not variance
     share specifically.

  3. Scenario/robustness analysis -- best-case / worst-case / baseline named
     parameter sets, evaluated as point simulations, to answer "how much do
     our qualitative conclusions change under a pessimistic vs optimistic
     reading of the evidence?"
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel


@dataclass
class PRCCResult:
    parameter_name: str
    prcc: float
    p_value: float


def partial_rank_correlation(sampled_params: dict[str, np.ndarray], output: np.ndarray) -> list[PRCCResult]:
    """Compute PRCC for each parameter in sampled_params against `output`,
    following the Marino et al. (2008) methodology: rank-transform all
    variables, then compute each parameter's partial correlation with the
    output controlling for all other sampled parameters.

    IMPORTANT -- pass only independently-sampled inputs, never a derived/
    composite quantity alongside its own constituent factors. E.g. if
    r_effective = r0 * contact_multiplier * intervention_multiplier, do NOT
    include r_effective in the same call as r0/contact_multiplier/
    intervention_multiplier: r_effective is then a deterministic (collinear)
    function of the other three, which breaks the partial-correlation
    residualization step and produces misleadingly small PRCCs for the true
    independent drivers (their effect gets silently absorbed into the
    composite variable instead). Diagnose this failure mode by checking
    whether a parameter's PRCC drops to ~0 with p >> 0.05 despite an
    obviously large one-way sensitivity swing (see one_way.py) for the same
    parameter -- that mismatch is the signature of exactly this issue, not
    genuine unimportance.
    """
    names = list(sampled_params.keys())
    if len(names) < 2:
        raise ValueError("PRCC requires at least 2 varying parameters to control for one another")

    # rank-transform everything (this is what makes it a *rank* correlation, robust to nonlinear-but-monotonic relationships)
    ranked = {name: stats.rankdata(arr) for name, arr in sampled_params.items()}
    ranked_output = stats.rankdata(output)

    results = []
    for target in names:
        others = [n for n in names if n != target]
        X_target = ranked[target]
        X_others = np.column_stack([ranked[n] for n in others]) if others else np.empty((len(X_target), 0))

        def residual(y: np.ndarray) -> np.ndarray:
            if X_others.shape[1] == 0:
                return y - y.mean()
            design = np.column_stack([np.ones(len(y)), X_others])
            coef, *_ = np.linalg.lstsq(design, y, rcond=None)
            return y - design @ coef

        resid_target = residual(X_target)
        resid_output = residual(ranked_output)
        r, p = stats.pearsonr(resid_target, resid_output)
        results.append(PRCCResult(parameter_name=target, prcc=float(r), p_value=float(p)))

    return sorted(results, key=lambda x: abs(x.prcc), reverse=True)


@dataclass
class LeaveOneOutResult:
    parameter_name: str
    full_variance: float
    variance_with_parameter_fixed: float
    variance_contribution_fraction: float


def leave_one_out_variance_contribution(
    population_size: int,
    initial_cases: int,
    baseline_point_estimates: dict[str, float],
    parameter_distributions: dict[str, "callable"],  # name -> scipy dist with .rvs(size, random_state)
    contact_multiplier: float,
    n_iterations: int = 2000,
    seed: int = 20260720,
) -> list[LeaveOneOutResult]:
    """For each parameter, compare output (attack rate) variance from a full
    Monte Carlo run (all parameters varying) against a run with that one
    parameter clamped to its point estimate -- the drop in variance
    attributes that share of total output uncertainty to the clamped
    parameter. This is the simulation-parameter analogue of "leave-one-out"
    (adapted from the meta-analysis leave-one-study-out convention referenced
    in this project's sensitivity-analysis-design skill guidance, applied
    here to simulation input parameters rather than included studies)."""

    def run_variance(fixed_param: str | None) -> float:
        rng = np.random.default_rng(seed)
        rates = np.empty(n_iterations)
        for i in range(n_iterations):
            draws = {}
            for name, dist in parameter_distributions.items():
                if name == fixed_param:
                    draws[name] = baseline_point_estimates[name]
                else:
                    draws[name] = float(dist.rvs(random_state=rng))
            r_eff = draws["r0"] * contact_multiplier
            cfg = BranchingProcessConfig(
                population_size=population_size, initial_cases=initial_cases,
                r_effective=r_eff, k_dispersion=draws["k_dispersion"],
            )
            rates[i] = BranchingProcessModel(cfg).run(rng).attack_rate
        return float(np.var(rates, ddof=1))

    full_var = run_variance(fixed_param=None)
    results = []
    for name in parameter_distributions:
        clamped_var = run_variance(fixed_param=name)
        contribution = (full_var - clamped_var) / full_var if full_var > 0 else 0.0
        results.append(LeaveOneOutResult(
            parameter_name=name, full_variance=full_var,
            variance_with_parameter_fixed=clamped_var,
            variance_contribution_fraction=contribution,
        ))
    return sorted(results, key=lambda r: r.variance_contribution_fraction, reverse=True)


@dataclass
class ScenarioComparisonResult:
    scenario_name: str
    parameters: dict[str, float]
    mean_attack_rate: float


def scenario_robustness_analysis(
    population_size: int,
    initial_cases: int,
    contact_multiplier: float,
    scenarios: dict[str, dict[str, float]],  # scenario_name -> {'r0':..., 'k_dispersion':...}
    n_reps: int = 1000,
    seed: int = 20260720,
) -> list[ScenarioComparisonResult]:
    """Best-case / worst-case / baseline (or any user-named set of) parameter
    combinations, each evaluated with n_reps stochastic realizations."""
    results = []
    for name, params in scenarios.items():
        rng = np.random.default_rng(seed)
        rates = np.empty(n_reps)
        for i in range(n_reps):
            cfg = BranchingProcessConfig(
                population_size=population_size, initial_cases=initial_cases,
                r_effective=params["r0"] * contact_multiplier, k_dispersion=params["k_dispersion"],
            )
            rates[i] = BranchingProcessModel(cfg).run(rng).attack_rate
        results.append(ScenarioComparisonResult(scenario_name=name, parameters=params, mean_attack_rate=float(rates.mean())))
    return results
