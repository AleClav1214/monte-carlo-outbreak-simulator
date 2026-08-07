"""
Data validation engine.

Implements requirement #12 (Data Validation) of the project spec:
  - schema validation        -> delegated to Pydantic models in schemas.py
  - parameter validation     -> range/domain checks, embedded in schemas.py validators
  - range checks             -> embedded in schemas.py validators
  - missing data handling    -> explicit policy below (fail loudly, never silently impute)
  - integrity checks         -> cross-file referential integrity, implemented here

Design principle: fail loudly, not silently.
A simulation run on a scenario with an undefined pathogen reference, or a
parameter set missing a required field, must raise before any Monte Carlo
iteration executes — not produce plausible-looking numbers from partially
wrong inputs. This mirrors the "validating-database-integrity" skill's
emphasis on referential integrity: never silently fill in what wasn't
provided.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from pydantic import ValidationError

from outbreak_simulator.data.schemas import OutbreakScenario, PathogenParameterSet


class DataIntegrityError(Exception):
    """Raised when loaded data fails schema, range, or referential-integrity checks."""


@dataclass
class ValidationReport:
    """Accumulates all problems found in a single validation pass, rather than
    failing on the first one — this makes fixing a bad data file far less
    tedious than a fail-fast-on-first-error design."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0

    def raise_if_invalid(self) -> None:
        if not self.is_valid:
            msg = "Data validation failed:\n" + "\n".join(f"  - {e}" for e in self.errors)
            raise DataIntegrityError(msg)

    def summary(self) -> str:
        lines = [f"Validation report: {len(self.errors)} error(s), {len(self.warnings)} warning(s)"]
        lines += [f"  ERROR: {e}" for e in self.errors]
        lines += [f"  WARNING: {w}" for w in self.warnings]
        return "\n".join(lines)


def validate_pathogen_dict(raw: dict, pathogen_id: str, report: ValidationReport) -> PathogenParameterSet | None:
    """Validate one pathogen's raw YAML dict against the schema."""
    try:
        return PathogenParameterSet.model_validate(raw)
    except ValidationError as exc:
        for err in exc.errors():
            loc = ".".join(str(x) for x in err["loc"])
            report.errors.append(f"pathogen '{pathogen_id}': {loc}: {err['msg']}")
        return None


def validate_scenario_dict(raw: dict, scenario_id: str, report: ValidationReport) -> OutbreakScenario | None:
    """Validate one scenario's raw YAML dict against the schema."""
    try:
        return OutbreakScenario.model_validate(raw)
    except ValidationError as exc:
        for err in exc.errors():
            loc = ".".join(str(x) for x in err["loc"])
            report.errors.append(f"scenario '{scenario_id}': {loc}: {err['msg']}")
        return None


def check_referential_integrity(
    pathogens: dict[str, PathogenParameterSet],
    scenarios: dict[str, OutbreakScenario],
    report: ValidationReport,
) -> None:
    """Cross-file checks that individual-file schema validation cannot catch."""
    for scenario_id, scenario in scenarios.items():
        if scenario.pathogen_id not in pathogens:
            report.errors.append(
                f"scenario '{scenario_id}' references unknown pathogen_id "
                f"'{scenario.pathogen_id}' (known: {sorted(pathogens.keys())})"
            )
        if not scenario.observed_outcomes:
            report.warnings.append(
                f"scenario '{scenario_id}' has no observed_outcomes — external validation "
                f"(requirement #9) will not be possible for this scenario until real-world "
                f"benchmark data is added"
            )
        if not scenario.limitations:
            report.warnings.append(
                f"scenario '{scenario_id}' declares no limitations — this is almost never "
                f"actually true and should be reviewed"
            )


def check_missing_data_policy(raw: dict, context: str, required_fields: set[str], report: ValidationReport) -> None:
    """Explicit missing-data policy: every required field must be *present*, even if the
    value is a documented placeholder — silent omission is what we are guarding against.
    This is stricter than Pydantic's own required-field check because it runs on the raw
    dict *before* schema coercion, so it also catches keys present with value `null` that
    a permissive schema might otherwise accept."""
    present = set(raw.keys())
    missing = required_fields - present
    if missing:
        report.errors.append(f"{context}: missing required fields (present but null, or absent): {missing}")
    for k in required_fields & present:
        if raw[k] is None:
            report.errors.append(f"{context}: field '{k}' is explicitly null — provide a value or remove the scenario")


def validate_all(pathogens_raw: dict[str, dict], scenarios_raw: dict[str, dict]) -> tuple[
    dict[str, PathogenParameterSet], dict[str, OutbreakScenario], ValidationReport
]:
    """Full validation pipeline: schema -> per-object -> referential integrity.

    Returns validated objects plus a report. Callers decide whether to call
    `report.raise_if_invalid()` immediately or inspect warnings first.
    """
    report = ValidationReport()

    pathogens: dict[str, PathogenParameterSet] = {}
    for pid, raw in pathogens_raw.items():
        required = {"pathogen_id", "parameters", "last_reviewed", "review_method"}
        check_missing_data_policy(raw, f"pathogen '{pid}'", required, report)
        validated = validate_pathogen_dict(raw, pid, report)
        if validated is not None:
            pathogens[pid] = validated

    scenarios: dict[str, OutbreakScenario] = {}
    for sid, raw in scenarios_raw.items():
        check_missing_data_policy(raw, f"scenario '{sid}'", {"scenario_id", "pathogen_id", "population"}, report)
        validated = validate_scenario_dict(raw, sid, report)
        if validated is not None:
            scenarios[sid] = validated

    check_referential_integrity(pathogens, scenarios, report)
    return pathogens, scenarios, report


def validate_files(parameters_path: Path, scenarios_path: Path) -> tuple[
    dict[str, PathogenParameterSet], dict[str, OutbreakScenario], ValidationReport
]:
    """Convenience wrapper: load YAML from disk, then run validate_all()."""
    import yaml  # local import keeps this module importable without PyYAML for pure schema use

    with open(parameters_path) as f:
        pathogens_raw = yaml.safe_load(f) or {}
    with open(scenarios_path) as f:
        scenarios_raw = yaml.safe_load(f) or {}
    return validate_all(pathogens_raw, scenarios_raw)
