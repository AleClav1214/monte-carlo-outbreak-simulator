"""
Scenario calibration workflow.

Applies the metrics in metrics.py to a scenario's actual Monte Carlo output
against its embedded real-world observed_outcomes (scenarios.yaml), and
reports results using the validation-tier framework documented in
docs/validation_plan.md: this module NEVER emits an unqualified "validated"
label. Every calibration result is tagged with which validation tier it
represents (internal resampling / external cohort / etc.) and an explicit
statement of what it does NOT establish.
"""

from __future__ import annotations

from dataclasses import dataclass

from outbreak_simulator.simulations.runner import ScenarioRunResult
from outbreak_simulator.validation.metrics import CoverageResult, posterior_predictive_check, predictive_coverage


@dataclass
class ObservationCalibration:
    observation_description: str
    observation_source: str
    coverage: CoverageResult
    ppc_pvalue: float
    interpretation: str


@dataclass
class ScenarioCalibrationReport:
    scenario_id: str
    validation_tier: str
    n_independent_benchmarks: int
    per_observation: list[ObservationCalibration]
    overall_statement: str
    NOT_established: list[str]


def calibrate_scenario(run_result: ScenarioRunResult) -> ScenarioCalibrationReport:
    """Compare a scenario's simulated attack-rate distribution against each
    of its embedded observed_outcomes."""
    scenario = run_result.scenario
    simulated = run_result.mc_result.raw_attack_rates
    per_obs = []

    for obs in scenario.observed_outcomes:
        if obs.attack_rate is None:
            continue
        cov = predictive_coverage(simulated, obs.attack_rate, interval_level=0.95)
        ppc = posterior_predictive_check(simulated, obs.attack_rate)

        if cov.covered and 10 <= cov.observed_percentile <= 90:
            interp = "Observed value is squarely within the model's central predictive mass."
        elif cov.covered:
            interp = (
                f"Observed value falls within the 95% predictive interval but in an extreme percentile "
                f"({cov.observed_percentile:.0f}th) -- consistent with the model, but only in its tail."
            )
        else:
            interp = (
                f"Observed value ({obs.attack_rate:.1%}) falls OUTSIDE the model's 95% predictive interval "
                f"{cov.predictive_interval} -- either this specific observation is an outlier even relative to "
                f"the modeled uncertainty, or a scenario parameter (contact_rate_multiplier, susceptible_fraction, "
                f"etc.) needs revision. Do not treat the model as calibrated to this benchmark without investigating."
            )
        per_obs.append(ObservationCalibration(
            observation_description=obs.description, observation_source=obs.source,
            coverage=cov, ppc_pvalue=ppc, interpretation=interp,
        ))

    n_benchmarks = len(per_obs)
    if n_benchmarks == 0:
        tier = "none (no external benchmark data available for this scenario)"
        overall = (
            "No external validation is possible for this scenario -- no observed_outcomes "
            "with a reported attack_rate are available."
        )
    elif n_benchmarks == 1:
        tier = "external cohort validation (n=1)"
        overall = (
            "A SINGLE real-world benchmark was checked. This can rule out gross miscalibration "
            "(the model producing outcomes wildly inconsistent with the one known real event) but cannot "
            "establish general external validity -- one data point cannot distinguish 'this setting is "
            "reliably modeled' from 'this specific case happened to be consistent with these parameters.'"
        )
    else:
        tier = f"external cohort validation (n={n_benchmarks}, still a small-sample check)"
        overall = (
            f"{n_benchmarks} independent real-world benchmarks were checked. This is a stronger check than "
            f"n=1 but still a small sample -- treat as suggestive consistency-checking, not formal external validation."
        )

    not_established = [
        "Functional validation: whether this model's outputs are reliable enough to inform real "
        "public-health decisions has NOT been assessed and is out of scope for this project.",
        "Translational/implementation validation: this tool has not been compared against, or reviewed "
        "alongside, established outbreak-modeling tools used operationally by public health agencies.",
        "Temporal validation: whether calibration holds across different time periods / pathogen "
        "variants / population-immunity states has not been assessed.",
    ]
    if n_benchmarks <= 1:
        not_established.insert(0, (
            "Site-based validation: calibration against multiple independent occurrences of "
            "this same scenario type has not been assessed (only one real-world instance is available)."
        ))

    return ScenarioCalibrationReport(
        scenario_id=scenario.scenario_id, validation_tier=tier, n_independent_benchmarks=n_benchmarks,
        per_observation=per_obs, overall_statement=overall, NOT_established=not_established,
    )


def print_calibration_report(report: ScenarioCalibrationReport) -> str:
    lines = [
        f"Calibration report: {report.scenario_id}",
        f"Validation tier: {report.validation_tier}",
        "",
        report.overall_statement,
        "",
    ]
    for obs in report.per_observation:
        lines.append(f"  Benchmark: {obs.observation_description}")
        lines.append(f"    Source: {obs.observation_source}")
        interval_str = tuple(round(x, 3) for x in obs.coverage.predictive_interval)
        lines.append(f"    Observed={obs.coverage.observed_value:.1%}, 95% predictive interval={interval_str}, "
                      f"percentile={obs.coverage.observed_percentile:.0f}, PPC p={obs.ppc_pvalue:.3f}")
        lines.append(f"    {obs.interpretation}")
        lines.append("")
    lines.append("Explicitly NOT established by this report:")
    for item in report.NOT_established:
        lines.append(f"  - {item}")
    return "\n".join(lines)
