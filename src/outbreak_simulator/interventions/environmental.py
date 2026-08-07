"""
Environmental interventions: masking, ventilation, air filtration.

Ventilation and air filtration are both modeled via the Wells-Riley
rebreathed-air-fraction mechanism (Wells 1955; standard reference
formulation: P_infection = 1 - exp(-I*q*p*t / Q), where Q is the
room's effective clean-air supply rate) -- both act on the same physical
quantity (effective air changes per hour, ACH), so are implemented through
one shared function, `ach_risk_multiplier()`, with air filtration treated as
equivalent additional ACH (a standard engineering equivalence: a HEPA unit
removing X m3/h of air at ~99.97% efficiency is treated as adding
approximately X/room_volume additional ACH -- see CDC's "clean air delivery
rate" framework).
"""

from __future__ import annotations

import numpy as np

from outbreak_simulator.interventions.base import Intervention


def ach_risk_multiplier(baseline_ach: float, improved_ach: float) -> float:
    """Wells-Riley-derived relative risk multiplier from increasing air
    changes per hour from baseline_ach to improved_ach, holding the quanta
    generation rate, breathing rate, and exposure time constant.

    Under Wells-Riley, steady-state infection risk for a given exposure
    time is approximately proportional to 1/Q (Q = outdoor/clean air supply
    rate, proportional to ACH for a fixed room volume) in the regime where
    the exponent I*q*p*t/Q is small (the common case for short/moderate
    exposures at reasonable ventilation) -- so risk_multiplier =
    baseline_ach / improved_ach approximates the relative risk reduction.
    This project uses the *approximation* rather than the full exponential
    form because scenario-level quanta generation rates (q) are not part of
    this project's evidence table (would require a dedicated aerosol-physics
    literature review) -- documented as a simplification, and flagged for
    sensitivity analysis when the exposure is long/high-risk (e.g. the choir
    scenario), where the linear approximation is weaker.
    """
    if baseline_ach <= 0 or improved_ach <= 0:
        raise ValueError("ACH values must be > 0")
    return float(np.clip(baseline_ach / improved_ach, 0.0, 1.0))


def masking(coverage: float = 0.7, evidence_basis: str = "observational") -> Intervention:
    """Mask-wearing intervention.

    Two named evidence bases are provided because the literature genuinely
    disagrees by effect size, not just by confidence interval width:
      - 'observational': pooled meta-analytic OR~0.35 (95% CI 0.24-0.51)
        for general/non-HCW populations, treating OR as an approximation to
        relative risk (reasonable when baseline attack rates are not
        extreme).
      - 'rct': RCT-only evidence is much more equivocal -- an adjusted OR
        around 0.92 (95% CI 0.81-1.04), not statistically significant, from
        pooled community-mask RCTs.
    Defaulting to 'observational' would overstate confidence if not paired
    with this caveat -- callers modeling a conservative/skeptical scenario
    should explicitly pass evidence_basis='rct'.
    """
    if evidence_basis == "observational":
        mult, low, high = 0.53, 0.36, 0.79  # non-HCW subgroup pooled OR, approximated as RR
        source = "Pooled observational-study meta-analysis, non-healthcare-worker subgroup (OR 0.53, 95% CI 0.36-0.79)"
    elif evidence_basis == "rct":
        mult, low, high = 0.92, 0.81, 1.04
        source = (
            "Pooled community-mask RCT evidence (adjusted OR 0.92, 95% CI 0.81-1.04) -- effect "
            "not statistically significant; included to represent genuine scientific "
            "disagreement, not to be treated as a confident estimate either"
        )
    else:
        raise ValueError("evidence_basis must be 'observational' or 'rct'")
    return Intervention(
        name=f"Masking ({evidence_basis} evidence base)",
        category="environmental",
        transmission_multiplier=mult,
        multiplier_low=low,
        multiplier_high=high,
        coverage=coverage,
        source=source,
        description=(
            "Per-contact transmission-probability reduction from mask-wearing; "
            "coverage = population compliance rate."
        ),
    )


def ventilation(baseline_ach: float, improved_ach: float, coverage: float = 1.0) -> Intervention:
    """Ventilation improvement (e.g. opening windows, upgrading HVAC).
    coverage=1.0 by default since improved room ventilation typically
    affects the whole shared airspace, not a compliance-dependent subset of
    individuals."""
    mult = ach_risk_multiplier(baseline_ach, improved_ach)
    return Intervention(
        name=f"Ventilation improvement ({baseline_ach:.1f} -> {improved_ach:.1f} ACH)",
        category="environmental",
        transmission_multiplier=mult,
        multiplier_low=max(mult * 0.7, 0.0),
        multiplier_high=min(mult * 1.3, 1.0),
        coverage=coverage,
        source=(
            "Wells-Riley rebreathed-air-fraction approximation (risk approx. proportional to "
            "1/ACH); CDC baseline (5 ACH) and healthcare airborne-isolation target "
            "(6-12 ACH) as reference points"
        ),
        description=(
            "Models ventilation as reducing the effective concentration of infectious "
            "aerosol via increased clean-air supply."
        ),
    )


def air_filtration(
    room_volume_m3: float, clean_air_delivery_rate_m3h: float, baseline_ach: float, coverage: float = 1.0
) -> Intervention:
    """Portable HEPA air filtration, modeled as equivalent additional ACH
    (CDC's Clean Air Delivery Rate framework: CADR / room_volume ~ additional ACH)."""
    if room_volume_m3 <= 0:
        raise ValueError("room_volume_m3 must be > 0")
    additional_ach = clean_air_delivery_rate_m3h / room_volume_m3
    improved_ach = baseline_ach + additional_ach
    mult = ach_risk_multiplier(baseline_ach, improved_ach)
    return Intervention(
        name=f"Portable HEPA filtration (+{additional_ach:.1f} equivalent ACH)",
        category="environmental",
        transmission_multiplier=mult,
        multiplier_low=max(mult * 0.7, 0.0),
        multiplier_high=min(mult * 1.3, 1.0),
        coverage=coverage,
        source=(
            "CDC Clean Air Delivery Rate (CADR) equivalence framework + Wells-Riley "
            "approximation (see ventilation() docstring)"
        ),
        description=(
            "Portable HEPA filtration treated as equivalent additional air changes per "
            "hour, stacked on top of baseline ventilation."
        ),
    )
