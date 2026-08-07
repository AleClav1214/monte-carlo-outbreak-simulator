"""
Intervention base abstraction.

Design principle: every intervention (vaccination, masking, ventilation,
testing, isolation, ...) is represented uniformly as a multiplicative
modifier on effective transmission risk, with an explicit `coverage` term
(the fraction of the population/contacts actually reached by the
intervention) separate from `efficacy` (how effective it is when it does
apply). This decomposition matters: a 90%-efficacious intervention with 20%
coverage has the same net effect as a 45%-efficacious one with 40% coverage
-- keeping these separate makes that visible and lets each be independently
uncertain (see intervention_multiplier_distribution below, used by
sensitivity/ and simulations/monte_carlo.py to propagate uncertainty in
BOTH efficacy and coverage, not just efficacy).

Composability: InterventionStack (stack.py) combines multiple simultaneous
interventions multiplicatively on the *hazard* (transmission-probability)
scale, which is the standard simplifying assumption for combining
independent partial protections (e.g. used in the Wells-Riley/dose-response
literature this project's ventilation model is based on) -- documented
explicitly as an assumption, not asserted as exact, since real interventions
can interact non-independently (e.g. masking matters less if ventilation
already dilutes aerosol to negligible concentration).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from scipy import stats


@dataclass
class Intervention:
    """A single intervention, represented as a transmission-risk modifier.

    transmission_multiplier: multiplies effective R (or beta) when the
        intervention is fully applied to a contact (1.0 = no effect,
        0.0 = completely blocks transmission). efficacy = 1 - multiplier.
    coverage: fraction of the population/contacts actually reached
        (compliance rate, vaccination coverage, fraction of air actually
        filtered, etc.) -- 1.0 = universal.
    multiplier_uncertainty_low/high: literature-reported uncertainty
        interval on transmission_multiplier, for Monte Carlo propagation.
    """

    name: str
    category: str  # "pharmaceutical" | "environmental" | "behavioral"
    transmission_multiplier: float
    multiplier_low: Optional[float] = None
    multiplier_high: Optional[float] = None
    coverage: float = 1.0
    source: str = ""
    description: str = ""

    def __post_init__(self):
        if not (0.0 <= self.transmission_multiplier <= 1.0):
            raise ValueError(f"{self.name}: transmission_multiplier must be in [0,1]")
        if not (0.0 <= self.coverage <= 1.0):
            raise ValueError(f"{self.name}: coverage must be in [0,1]")

    @property
    def efficacy(self) -> float:
        return 1.0 - self.transmission_multiplier

    def net_multiplier(self) -> float:
        """Population-level multiplier accounting for partial coverage:
        uncovered contacts are unaffected, covered contacts get the full
        transmission_multiplier -- net = 1 - coverage*(1 - multiplier)."""
        return 1.0 - self.coverage * (1.0 - self.transmission_multiplier)

    def multiplier_distribution(self):
        """A Beta-ish sampling distribution over transmission_multiplier for
        Monte Carlo propagation of intervention-effectiveness uncertainty,
        built the same way models/distributions.py builds parameter
        distributions from an evidence-table entry."""
        from outbreak_simulator.models.distributions import beta_from_mean_ci

        if self.multiplier_low is None or self.multiplier_high is None:
            return stats.uniform(loc=self.transmission_multiplier, scale=0.0)  # degenerate: no reported uncertainty
        mean = max(min(self.transmission_multiplier, 0.999), 0.001)
        return beta_from_mean_ci(mean, max(self.multiplier_low, 1e-6), min(self.multiplier_high, 1 - 1e-6))
