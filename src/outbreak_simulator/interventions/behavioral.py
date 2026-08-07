"""
Behavioral/administrative interventions: testing, isolation, quarantine,
occupancy reduction.

These interventions act on transmission risk through detection-and-removal
(testing + isolation) or contact-reduction (quarantine of contacts,
occupancy limits) mechanisms rather than per-contact transmission
probability -- but are represented through the same Intervention interface
so they compose identically with the masking/ventilation/vaccination
modifiers in InterventionStack.
"""

from __future__ import annotations

from outbreak_simulator.interventions.base import Intervention


def testing_isolation(
    test_sensitivity: float,
    testing_frequency_per_infectious_period: float,
    isolation_effectiveness: float = 0.9,
    coverage: float = 1.0,
) -> Intervention:
    """Routine/surveillance testing that leads to isolation of detected cases.

    Approximates the probability a case is detected and isolated before
    completing its full transmission potential as
    1 - (1 - test_sensitivity)^testing_frequency_per_infectious_period,
    then scales transmission by isolation_effectiveness once detected
    (isolation is rarely 100% effective at preventing all further
    transmission -- e.g. household contacts, essential activities).

    Parameters
    ----------
    test_sensitivity: per-test probability of correctly detecting an
        infectious case (varies hugely by test type/timing -- pass the
        value appropriate to the assay in use; this project does not
        hardcode a default since it is too context-dependent to be a
        reasonable one-size-fits-all assumption).
    testing_frequency_per_infectious_period: expected number of tests
        administered during a case's infectious period (e.g. daily testing
        over a 5-day infectious period ~= 5).
    isolation_effectiveness: fraction of a detected case's remaining
        transmission potential that isolation actually prevents.
    """
    if not (0 <= test_sensitivity <= 1):
        raise ValueError("test_sensitivity must be in [0,1]")
    if testing_frequency_per_infectious_period < 0:
        raise ValueError("testing_frequency_per_infectious_period must be >= 0")
    detection_probability = 1 - (1 - test_sensitivity) ** testing_frequency_per_infectious_period
    net_reduction = detection_probability * isolation_effectiveness
    mult = 1 - net_reduction
    return Intervention(
        name=f"Testing + isolation (sensitivity={test_sensitivity:.2f}, freq={testing_frequency_per_infectious_period:.1f}/infectious period)",
        category="behavioral",
        transmission_multiplier=max(mult, 0.0),
        multiplier_low=max(mult - 0.1, 0.0),
        multiplier_high=min(mult + 0.1, 1.0),
        coverage=coverage,
        source="Standard test-and-isolate epidemiological model: P(detect) = 1-(1-sensitivity)^n_tests, scaled by isolation effectiveness",
        description="Reduces transmission by detecting and isolating infectious individuals before they complete their transmission potential.",
    )


def quarantine_contacts(quarantine_effectiveness: float, contact_tracing_coverage: float) -> Intervention:
    """Quarantine of TRACED CONTACTS of a known case (distinct from isolation
    of the case itself). quarantine_effectiveness = fraction reduction in
    onward transmission achieved for contacts who are successfully
    quarantined (accounts for imperfect adherence and the possibility of
    transmission before quarantine begins); contact_tracing_coverage =
    fraction of true contacts who are actually identified and quarantined.
    """
    if not (0 <= quarantine_effectiveness <= 1):
        raise ValueError("quarantine_effectiveness must be in [0,1]")
    return Intervention(
        name="Contact quarantine",
        category="behavioral",
        transmission_multiplier=1 - quarantine_effectiveness,
        coverage=contact_tracing_coverage,
        source="Standard contact-tracing-and-quarantine model: effect requires both successful tracing (coverage) and adherent quarantine (effectiveness)",
        description="Reduces onward transmission from contacts of known cases who are identified and quarantined before becoming infectious.",
    )


def occupancy_reduction(occupancy_fraction: float) -> Intervention:
    """Reducing the number of people present (e.g. capping event attendance,
    splitting a cohort into smaller groups, hybrid/remote arrangements).

    Modeled as an approximately linear reduction in effective transmission
    for a well-mixed setting: halving occupancy roughly halves the number
    of close-contact transmission opportunities per unit time. This is a
    simplification (real contact structure is rarely perfectly well-mixed,
    and some settings have a fixed number of unavoidable close contacts
    regardless of total occupancy -- e.g. roommates) -- flagged as such.
    """
    if not (0 < occupancy_fraction <= 1):
        raise ValueError("occupancy_fraction must be in (0,1]")
    return Intervention(
        name=f"Occupancy reduction (to {occupancy_fraction:.0%} of baseline)",
        category="behavioral",
        transmission_multiplier=occupancy_fraction,
        multiplier_low=max(occupancy_fraction - 0.15, 0.0),
        multiplier_high=min(occupancy_fraction + 0.15, 1.0),
        coverage=1.0,
        source="Well-mixed-population approximation: contact opportunities per unit time scale approximately linearly with occupancy",
        description="Reduces transmission by reducing the number of people sharing the space, approximated as linear in a well-mixed setting.",
    )
