"""
Intervention stacking and scenario comparison.

Combines multiple simultaneous interventions on the hazard (transmission
probability) scale -- the standard, explicitly-flagged simplifying
assumption for combining independent partial protections (see
interventions/base.py module docstring for the caveat about non-independent
interactions).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from outbreak_simulator.interventions.base import Intervention


@dataclass
class InterventionStack:
    """A named, ordered collection of simultaneously-applied interventions."""

    name: str
    interventions: list[Intervention]

    def combined_multiplier(self) -> float:
        """Multiplicative combination of each intervention's net (coverage-
        adjusted) multiplier onto baseline transmission risk."""
        m = 1.0
        for iv in self.interventions:
            m *= iv.net_multiplier()
        return m

    def sample_combined_multiplier(self, rng: np.random.Generator) -> float:
        """Monte Carlo draw of the combined multiplier, propagating
        uncertainty in EACH intervention's effectiveness (not just using
        each one's point estimate) -- coverage is treated as fixed/known
        (a policy choice, not usually subject to the same kind of
        measurement uncertainty as an effect size) while transmission_multiplier
        is resampled from its literature-reported interval each iteration."""
        m = 1.0
        for iv in self.interventions:
            sampled_mult = float(iv.multiplier_distribution().rvs(random_state=rng))
            sampled_mult = float(np.clip(sampled_mult, 0.0, 1.0))
            net = 1.0 - iv.coverage * (1.0 - sampled_mult)
            m *= net
        return m

    def summary_table(self) -> list[dict]:
        rows = []
        for iv in self.interventions:
            rows.append({
                "name": iv.name,
                "category": iv.category,
                "efficacy": iv.efficacy,
                "coverage": iv.coverage,
                "net_multiplier": iv.net_multiplier(),
                "source": iv.source,
            })
        return rows


def compare_scenarios(baseline_r_effective: float, stacks: dict[str, InterventionStack]) -> dict[str, float]:
    """Support for requirement #6's 'scenario comparison': given a baseline
    effective R and a set of named intervention stacks (including,
    conventionally, a 'no_intervention' stack with an empty interventions
    list), return the resulting effective R under each."""
    result = {}
    for name, stack in stacks.items():
        result[name] = baseline_r_effective * stack.combined_multiplier()
    return result


def no_intervention() -> InterventionStack:
    return InterventionStack(name="No intervention (baseline)", interventions=[])
