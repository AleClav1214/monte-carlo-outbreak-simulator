"""
Pharmaceutical interventions: vaccination.

Vaccine effectiveness against ONWARD TRANSMISSION (what this module models)
is a distinct, generally smaller and much more heterogeneous quantity than
vaccine effectiveness against acquiring infection or against severe disease
-- conflating these is a common error this project deliberately avoids. See
each factory function's `source` for the specific literature this project
draws from.
"""

from __future__ import annotations

from outbreak_simulator.interventions.base import Intervention


def measles_vaccination(coverage: float = 0.95) -> Intervention:
    """Single-dose MMR, using the 1974 Sioux Reservation outbreak-derived
    efficacy estimate (97.3%, 95% CI 80.1-99.9%) -- see
    data/parameters/pathogens.yaml: measles.vaccine_efficacy.
    `coverage` defaults to 0.95, the commonly cited herd-immunity threshold
    for measles given R0~12-18, but should be set to the scenario's actual
    assumed/observed coverage."""
    return Intervention(
        name="MMR vaccination (single dose, historical estimate)",
        category="pharmaceutical",
        transmission_multiplier=1 - 0.973,
        multiplier_low=1 - 0.999,
        multiplier_high=1 - 0.801,
        coverage=coverage,
        source="Cherry et al., 1974 Sioux Reservation measles epidemic (see pathogens.yaml)",
        description="Reduces susceptibility to infection (and thus onward transmission) among vaccinated individuals.",
    )


def sars_cov_2_vaccination(coverage: float = 0.7, variant_era: str = "delta") -> Intervention:
    """COVID-19 vaccine effectiveness against ONWARD transmission (not
    against acquiring infection, which is a different, generally higher
    figure). Delta-era: ~63% reduction in transmission to unvaccinated
    contacts (Netherlands household study, Eyre et al.-adjacent 2022
    Delta-era estimates). Effectiveness against transmission was
    substantially lower and more variable for Omicron -- pass
    variant_era='omicron' for that regime."""
    if variant_era == "delta":
        mult, low, high = 1 - 0.50, 1 - 0.63, 1 - 0.20  # central ~50%, bounded by the 20-63% literature range found
        source = "Netherlands household-contact study, Delta era (63% reduction in transmission to unvaccinated contacts; range reflects broader 20-63% literature)"
    elif variant_era == "omicron":
        mult, low, high = 1 - 0.35, 1 - 0.60, 1 - 0.10
        source = "Multiple Omicron-era studies reporting substantially lower and more heterogeneous transmission-blocking effectiveness (range ~10-60% depending on booster timing and age)"
    else:
        raise ValueError("variant_era must be 'delta' or 'omicron'")
    return Intervention(
        name=f"SARS-CoV-2 vaccination ({variant_era} era, effectiveness against onward transmission)",
        category="pharmaceutical",
        transmission_multiplier=mult,
        multiplier_low=low,
        multiplier_high=high,
        coverage=coverage,
        source=source,
        description="Effectiveness specifically against transmitting to others if infected -- distinct from (and generally lower than) effectiveness against acquiring infection.",
    )


def generic_vaccination(efficacy: float, efficacy_low: float, efficacy_high: float, coverage: float, source: str) -> Intervention:
    """Construct a vaccination intervention for a pathogen without a
    dedicated factory function above -- requires the caller to supply a
    literature-sourced efficacy estimate explicitly rather than defaulting
    to an unsourced guess."""
    return Intervention(
        name="Vaccination (user-specified)",
        category="pharmaceutical",
        transmission_multiplier=1 - efficacy,
        multiplier_low=1 - efficacy_high,
        multiplier_high=1 - efficacy_low,
        coverage=coverage,
        source=source,
        description="User-supplied vaccine effectiveness against onward transmission.",
    )
