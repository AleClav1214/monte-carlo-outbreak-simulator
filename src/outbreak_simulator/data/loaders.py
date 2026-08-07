"""
Convenience loaders for pathogen parameters and scenarios.

These wrap data.validation.validate_files() with sensible defaults (the
package's bundled data files) and small ergonomic helpers used throughout
models/, simulations/, sensitivity/, and visualization/.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from outbreak_simulator.data.schemas import OutbreakScenario, ParameterEstimate, PathogenParameterSet
from outbreak_simulator.data.validation import validate_files

_PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_PARAMETERS_PATH = _PACKAGE_ROOT / "parameters" / "pathogens.yaml"
DEFAULT_SCENARIOS_PATH = _PACKAGE_ROOT / "scenarios" / "scenarios.yaml"
# NOTE: real-world benchmark outcomes are embedded directly on each scenario
# (OutbreakScenario.observed_outcomes) rather than kept in a separate file --
# this avoids a second file that could silently drift out of sync with the
# scenario it describes. See scenarios.yaml.


@lru_cache(maxsize=4)
def load_all(
    parameters_path: Path = DEFAULT_PARAMETERS_PATH,
    scenarios_path: Path = DEFAULT_SCENARIOS_PATH,
    strict: bool = True,
) -> tuple[dict[str, PathogenParameterSet], dict[str, OutbreakScenario]]:
    """Load and validate the bundled (or custom) pathogen + scenario data.

    Parameters
    ----------
    strict:
        If True (default), raise DataIntegrityError on any validation error.
        If False, return whatever validated successfully and print warnings
        (useful in exploratory notebooks; never use strict=False in
        production pipelines or CI).
    """
    pathogens, scenarios, report = validate_files(parameters_path, scenarios_path)
    if strict:
        report.raise_if_invalid()
    elif not report.is_valid:
        print(report.summary())
    if report.warnings:
        for w in report.warnings:
            print(f"[data warning] {w}")
    return pathogens, scenarios


def get_pathogen(pathogen_id: str) -> PathogenParameterSet:
    pathogens, _ = load_all()
    if pathogen_id not in pathogens:
        raise KeyError(f"Unknown pathogen_id '{pathogen_id}'. Known: {sorted(pathogens.keys())}")
    return pathogens[pathogen_id]


def get_scenario(scenario_id: str) -> OutbreakScenario:
    _, scenarios = load_all()
    if scenario_id not in scenarios:
        raise KeyError(f"Unknown scenario_id '{scenario_id}'. Known: {sorted(scenarios.keys())}")
    return scenarios[scenario_id]


def get_parameter(pathogen_id: str, parameter_name: str) -> ParameterEstimate:
    pathogen = get_pathogen(pathogen_id)
    if parameter_name not in pathogen.parameters:
        raise KeyError(
            f"Pathogen '{pathogen_id}' has no parameter '{parameter_name}'. "
            f"Known: {sorted(pathogen.parameters.keys())}"
        )
    return pathogen.parameters[parameter_name]


def list_pathogens() -> list[str]:
    pathogens, _ = load_all()
    return sorted(pathogens.keys())


def list_scenarios() -> list[str]:
    _, scenarios = load_all()
    return sorted(scenarios.keys())
