"""Unit tests: data schemas and validation engine."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from outbreak_simulator.data.loaders import get_parameter, get_pathogen, get_scenario, list_pathogens, list_scenarios
from outbreak_simulator.data.schemas import DistributionFamily, EvidenceQuality, ParameterEstimate
from outbreak_simulator.data.validation import (
    DataIntegrityError,
    ValidationReport,
    validate_all,
)


class TestParameterEstimateSchema:
    def _valid_kwargs(self, **overrides):
        base = dict(
            name="r0", display_name="R0", distribution=DistributionFamily.GAMMA,
            point_estimate=2.5, low=2.0, high=3.0, unit="dimensionless",
            source="Some Study 2020", evidence_quality=EvidenceQuality.SINGLE_STUDY,
            justification="Because the study said so and this justification is long enough.",
        )
        base.update(overrides)
        return base

    def test_valid_parameter_constructs(self):
        p = ParameterEstimate(**self._valid_kwargs())
        assert p.point_estimate == 2.5

    def test_rejects_point_estimate_outside_interval(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(point_estimate=10.0))

    def test_rejects_low_greater_than_high(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(low=5.0, high=1.0, point_estimate=3.0))

    def test_rejects_out_of_range_probability(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(unit="probability", point_estimate=1.5, low=1.0, high=2.0))

    def test_rejects_negative_duration(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(unit="days", point_estimate=-1.0, low=-2.0, high=0.0))

    def test_rejects_short_justification(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(justification="too short"))

    def test_rejects_missing_source(self):
        with pytest.raises(ValidationError):
            ParameterEstimate(**self._valid_kwargs(source=""))


class TestBundledData:
    def test_all_pathogens_load(self):
        pathogens = list_pathogens()
        assert set(pathogens) == {"sars_cov_2", "mpox", "influenza", "norovirus", "measles", "varicella"}

    def test_all_scenarios_load(self):
        scenarios = list_scenarios()
        assert len(scenarios) == 8

    def test_every_pathogen_has_required_core_parameters(self):
        required = {"r0", "incubation_period", "infectious_period", "k_dispersion", "secondary_attack_rate"}
        for pid in list_pathogens():
            pathogen = get_pathogen(pid)
            assert required.issubset(set(pathogen.parameters.keys())), f"{pid} missing required parameters"

    def test_every_scenario_references_valid_pathogen(self):
        for sid in list_scenarios():
            scenario = get_scenario(sid)
            get_pathogen(scenario.pathogen_id)  # raises KeyError if invalid -- test fails if it does

    def test_get_parameter_convenience_function(self):
        r0 = get_parameter("sars_cov_2", "r0")
        assert 2.0 < r0.point_estimate < 3.5

    def test_unknown_pathogen_raises_keyerror(self):
        with pytest.raises(KeyError):
            get_pathogen("nonexistent_pathogen")

    def test_low_confidence_parameters_are_flagged(self):
        """Regression test: this project explicitly flags weaker-evidence
        parameters (e.g. norovirus R0, several k_dispersion placeholders)
        rather than presenting them with the same confidence as
        meta-analytic figures -- verify at least one such flag exists and
        that the mechanism itself works."""
        norovirus = get_pathogen("norovirus")
        assert norovirus.parameters["r0"].evidence_quality == EvidenceQuality.LOW_CONFIDENCE


class TestReferentialIntegrity:
    def test_catches_scenario_with_unknown_pathogen(self):
        pathogens_raw = {"flu": {
            "pathogen_id": "flu", "display_name": "Flu", "pathogen_class": "respiratory",
            "last_reviewed": "2026-01-01", "review_method": "test",
            "parameters": {
                name: {
                    "name": name, "display_name": name, "distribution": "gamma",
                    "point_estimate": 1.0, "low": 0.5, "high": 1.5, "unit": "dimensionless",
                    "source": "test source", "evidence_quality": "single_study",
                    "justification": "test justification long enough to pass validation",
                } for name in ["r0", "incubation_period", "infectious_period", "k_dispersion", "secondary_attack_rate"]
            },
        }}
        scenarios_raw = {"bad_scenario": {
            "scenario_id": "bad_scenario", "display_name": "Bad", "pathogen_id": "does_not_exist",
            "description": "test", "population": {"population_size": 10, "setting_type": "test"},
        }}
        pathogens, scenarios, report = validate_all(pathogens_raw, scenarios_raw)
        assert not report.is_valid
        assert any("unknown pathogen_id" in e for e in report.errors)

    def test_report_raises_on_invalid(self):
        report = ValidationReport(errors=["fake error"])
        with pytest.raises(DataIntegrityError):
            report.raise_if_invalid()

    def test_report_does_not_raise_when_valid(self):
        report = ValidationReport()
        report.raise_if_invalid()  # should not raise
