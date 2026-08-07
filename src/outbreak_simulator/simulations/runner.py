"""
High-level scenario runner.

This is the main user-facing entry point: given a scenario_id (and
optionally an InterventionStack), it handles loading validated data,
constructing the appropriate model (branching process, by default, since
all bundled scenarios are closed-setting/superspreading-relevant; SEIR is
available for larger/longer-duration analyses), composing intervention
effects, and running the full Monte Carlo batch -- see examples/ for
end-to-end usage.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from outbreak_simulator.data import OutbreakScenario, PathogenParameterSet, get_pathogen, get_scenario
from outbreak_simulator.interventions.stack import InterventionStack, no_intervention
from outbreak_simulator.models.branching_process import BranchingProcessConfig, BranchingProcessModel
from outbreak_simulator.models.distributions import sample_parameter
from outbreak_simulator.simulations.convergence import ConvergenceReport, assess_convergence
from outbreak_simulator.simulations.monte_carlo import MonteCarloConfig, MonteCarloResult, run_monte_carlo


@dataclass
class ScenarioRunResult:
    scenario: OutbreakScenario
    pathogen: PathogenParameterSet
    intervention_stack: InterventionStack
    mc_result: MonteCarloResult
    convergence: ConvergenceReport


def run_scenario(
    scenario_id: str,
    intervention_stack: InterventionStack | None = None,
    n_iterations: int = 5000,
    seed: int = 20260720,
    susceptible_fraction: float = 1.0,
    store_results: bool = False,
) -> ScenarioRunResult:
    """Run a bundled scenario through the branching process model with full
    Monte Carlo uncertainty propagation over (r0, k_dispersion, and -- if
    provided -- intervention effectiveness).

    susceptible_fraction: fraction of the scenario's population_size assumed
    susceptible at baseline (1.0 = fully susceptible). This matters
    especially for the measles/varicella school scenarios, where population
    vaccination coverage determines how much of the nominal population_size
    is actually at risk -- see scenarios.yaml assumptions for those
    scenarios, which deliberately do not hardcode this.
    """
    scenario = get_scenario(scenario_id)
    pathogen = get_pathogen(scenario.pathogen_id)
    stack = intervention_stack if intervention_stack is not None else no_intervention()

    r0_param = pathogen.parameters["r0"]
    k_param = pathogen.parameters["k_dispersion"]
    effective_population = max(int(round(scenario.population.population_size * susceptible_fraction)), scenario.population.initial_cases)

    def factory(rng: np.random.Generator):
        r0 = float(sample_parameter(r0_param, rng))
        k = float(sample_parameter(k_param, rng))
        contact_mult = scenario.population.contact_rate_multiplier
        if scenario.population.contact_rate_multiplier_low is not None and scenario.population.contact_rate_multiplier_high is not None:
            # propagate uncertainty in the setting's own contact_rate_multiplier, not just the pathogen's r0/k
            low, high = scenario.population.contact_rate_multiplier_low, scenario.population.contact_rate_multiplier_high
            contact_mult = float(rng.triangular(low, scenario.population.contact_rate_multiplier, high))
        intervention_mult = stack.sample_combined_multiplier(rng)
        r_effective = r0 * contact_mult * intervention_mult

        cfg = BranchingProcessConfig(
            population_size=effective_population,
            initial_cases=min(scenario.population.initial_cases, effective_population),
            r_effective=r_effective,
            k_dispersion=k,
        )
        model = BranchingProcessModel(cfg)
        return model, {
            "r0": r0, "k_dispersion": k, "contact_rate_multiplier": contact_mult,
            "intervention_multiplier": intervention_mult, "r_effective": r_effective,
        }

    mc_config = MonteCarloConfig(n_iterations=n_iterations, master_seed=seed)
    mc_result = run_monte_carlo(factory, mc_config, store_results=store_results)
    convergence = assess_convergence(mc_result.raw_attack_rates)

    return ScenarioRunResult(
        scenario=scenario, pathogen=pathogen, intervention_stack=stack,
        mc_result=mc_result, convergence=convergence,
    )


def print_summary(result: ScenarioRunResult) -> str:
    """Human-readable summary of a scenario run -- used by examples/ and CLI-style usage."""
    s = result.mc_result.attack_rate_summary
    lines = [
        f"Scenario: {result.scenario.display_name}  (pathogen: {result.pathogen.display_name})",
        f"Intervention stack: {result.intervention_stack.name}",
        f"Monte Carlo iterations: {result.mc_result.config.n_iterations} (seed={result.mc_result.config.master_seed})",
        "",
        f"Attack rate  -- mean {s.mean:.1%}, median {s.median:.1%}, 95% UI [{s.percentile_5:.1%}, {s.percentile_95:.1%}]",
        f"Extinction probability: {result.mc_result.extinction_probability:.1%}",
        f"Convergence: {result.convergence.recommendation}",
    ]
    if result.scenario.observed_outcomes:
        lines.append("")
        lines.append("Observed real-world benchmark(s):")
        for obs in result.scenario.observed_outcomes:
            if obs.attack_rate is not None:
                lines.append(f"  - {obs.description}: attack rate {obs.attack_rate:.1%} ({obs.source})")
    return "\n".join(lines)
