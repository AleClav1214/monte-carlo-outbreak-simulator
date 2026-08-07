"""Unit tests: interventions/ (individual modifiers, stacking, scenario comparison)."""

from __future__ import annotations

import numpy as np
import pytest

from outbreak_simulator.interventions.base import Intervention
from outbreak_simulator.interventions.behavioral import occupancy_reduction, quarantine_contacts
from outbreak_simulator.interventions.behavioral import testing_isolation as apply_testing_isolation
from outbreak_simulator.interventions.environmental import ach_risk_multiplier, masking, ventilation
from outbreak_simulator.interventions.pharmaceutical import measles_vaccination, sars_cov_2_vaccination
from outbreak_simulator.interventions.stack import InterventionStack, compare_scenarios, no_intervention


class TestInterventionBase:
    def test_rejects_multiplier_outside_unit_interval(self):
        with pytest.raises(ValueError):
            Intervention(name="bad", category="test", transmission_multiplier=1.5)

    def test_net_multiplier_full_coverage_equals_transmission_multiplier(self):
        iv = Intervention(name="test", category="test", transmission_multiplier=0.5, coverage=1.0)
        assert iv.net_multiplier() == 0.5

    def test_net_multiplier_zero_coverage_equals_one(self):
        iv = Intervention(name="test", category="test", transmission_multiplier=0.1, coverage=0.0)
        assert iv.net_multiplier() == 1.0  # no effect if nobody is covered

    def test_partial_coverage_dilutes_effect(self):
        iv = Intervention(name="test", category="test", transmission_multiplier=0.5, coverage=0.5)
        assert iv.net_multiplier() == pytest.approx(0.75)  # 1 - 0.5*(1-0.5)


class TestEnvironmentalInterventions:
    def test_ach_risk_multiplier_higher_ach_lower_risk(self):
        assert ach_risk_multiplier(1.0, 6.0) < ach_risk_multiplier(1.0, 2.0)

    def test_ach_risk_multiplier_equal_ach_no_effect(self):
        assert ach_risk_multiplier(4.0, 4.0) == pytest.approx(1.0)

    def test_ach_rejects_nonpositive(self):
        with pytest.raises(ValueError):
            ach_risk_multiplier(0.0, 5.0)

    def test_masking_observational_stronger_than_rct(self):
        obs = masking(coverage=1.0, evidence_basis="observational")
        rct = masking(coverage=1.0, evidence_basis="rct")
        assert obs.transmission_multiplier < rct.transmission_multiplier

    def test_masking_rejects_unknown_evidence_basis(self):
        with pytest.raises(ValueError):
            masking(evidence_basis="not_a_real_option")

    def test_ventilation_intervention_has_expected_multiplier(self):
        iv = ventilation(baseline_ach=2.0, improved_ach=8.0)
        assert iv.transmission_multiplier == pytest.approx(0.25)


class TestPharmaceuticalInterventions:
    def test_measles_vaccination_high_efficacy(self):
        iv = measles_vaccination(coverage=0.95)
        assert iv.efficacy > 0.9

    def test_covid_vaccination_omicron_weaker_than_delta(self):
        delta = sars_cov_2_vaccination(coverage=1.0, variant_era="delta")
        omicron = sars_cov_2_vaccination(coverage=1.0, variant_era="omicron")
        assert omicron.efficacy < delta.efficacy

    def test_covid_vaccination_rejects_unknown_variant(self):
        with pytest.raises(ValueError):
            sars_cov_2_vaccination(variant_era="not_a_variant")


class TestBehavioralInterventions:
    def test_testing_isolation_more_frequent_testing_more_effective(self):
        low_freq = apply_testing_isolation(test_sensitivity=0.8, testing_frequency_per_infectious_period=1)
        high_freq = apply_testing_isolation(test_sensitivity=0.8, testing_frequency_per_infectious_period=5)
        assert high_freq.transmission_multiplier < low_freq.transmission_multiplier

    def test_occupancy_reduction_halving_halves_multiplier(self):
        iv = occupancy_reduction(occupancy_fraction=0.5)
        assert iv.transmission_multiplier == pytest.approx(0.5)

    def test_occupancy_reduction_rejects_out_of_range(self):
        with pytest.raises(ValueError):
            occupancy_reduction(occupancy_fraction=1.5)

    def test_quarantine_requires_both_tracing_and_effectiveness(self):
        iv = quarantine_contacts(quarantine_effectiveness=0.8, contact_tracing_coverage=0.3)
        # low tracing coverage should limit net effect even with high per-quarantine effectiveness
        assert iv.net_multiplier() > 0.7


class TestInterventionStack:
    def test_empty_stack_no_effect(self):
        stack = no_intervention()
        assert stack.combined_multiplier() == pytest.approx(1.0)

    def test_stack_multiplies_effects(self):
        iv1 = Intervention(name="a", category="test", transmission_multiplier=0.5, coverage=1.0)
        iv2 = Intervention(name="b", category="test", transmission_multiplier=0.5, coverage=1.0)
        stack = InterventionStack("test", [iv1, iv2])
        assert stack.combined_multiplier() == pytest.approx(0.25)

    def test_compare_scenarios_orders_correctly(self):
        strong = InterventionStack("strong", [Intervention(name="a", category="t", transmission_multiplier=0.1)])
        weak = InterventionStack("weak", [Intervention(name="a", category="t", transmission_multiplier=0.8)])
        result = compare_scenarios(5.0, {"baseline": no_intervention(), "strong": strong, "weak": weak})
        assert result["strong"] < result["weak"] < result["baseline"]

    def test_sample_combined_multiplier_reproducible(self):
        iv = Intervention(name="a", category="t", transmission_multiplier=0.5, multiplier_low=0.3, multiplier_high=0.7)
        stack = InterventionStack("test", [iv])
        v1 = stack.sample_combined_multiplier(np.random.default_rng(7))
        v2 = stack.sample_combined_multiplier(np.random.default_rng(7))
        assert v1 == v2
