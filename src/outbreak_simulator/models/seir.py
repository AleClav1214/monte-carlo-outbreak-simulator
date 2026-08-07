"""
Stochastic SEIR model, implemented as a generic reaction-network engine.

Architectural decision: rather than hardcoding an S/E/I/R state vector and
three fixed transitions, this module implements a general stochastic
reaction-network simulator (Gillespie's Direct Method / SSA) parameterized
by an arbitrary list of Reactions, each with a rate function and a state
delta. `build_seir_reactions()` below is the default SEIR instantiation of
this engine, but the same engine directly supports the project's
"configurable compartments" requirement -- adding an asymptomatic split,
a quarantine compartment, or a hospitalized compartment is a matter of
adding Reactions, not writing a new simulation loop.

Why Gillespie SSA rather than tau-leaping: the outbreak settings this
project models (choirs, barracks, dormitories, schools) have populations in
the tens to low thousands -- small enough that exact SSA is computationally
fine, and SSA avoids the tau-leaping approximation error that mainly exists
to buy speed on much larger populations than these scenarios need. This is
a deliberate scale-appropriate choice, documented as such in
docs/architecture.md.

Intervention hooks: interventions are represented as time-varying or
state-dependent multipliers on reaction rates (see interventions/base.py
for the Intervention protocol) rather than as special-cased branches in the
simulation loop -- an intervention is just a function of (state, t) that
returns a rate multiplier, applied uniformly to whichever reaction(s) it
targets.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from outbreak_simulator.models.base import OutbreakModel, SimulationResult

# A reaction's rate_fn takes (state: dict[str, int], t: float) -> float (a non-negative rate).
# Its delta is the state change applied when the reaction fires: dict[str, int].
RateFn = Callable[[dict, float], float]


@dataclass
class Reaction:
    name: str
    rate_fn: RateFn
    delta: dict[str, int]


@dataclass
class SEIRConfig:
    population_size: int
    initial_exposed: int = 0
    initial_infectious: int = 1
    beta: float = 0.5  # transmission rate (contacts/time * probability of transmission per contact)
    sigma: float = 0.2  # E -> I progression rate (1 / mean incubation-ish latent period)
    gamma: float = 0.14  # I -> R recovery rate (1 / mean infectious period)
    max_time: float = 365.0
    max_events: int = 2_000_000  # safety cap against runaway loops
    intervention_beta_multiplier: Callable[[dict, float], float] | None = None  # (state, t) -> multiplier on beta
    extra_reactions: list[Reaction] = field(default_factory=list)  # for configurable-compartment extensions


def build_seir_reactions(cfg: SEIRConfig) -> list[Reaction]:
    """Default SEIR reaction set: S->E (infection), E->I (progression), I->R (recovery)."""

    def infection_rate(state: dict, t: float) -> float:
        n = cfg.population_size
        beta_t = cfg.beta
        if cfg.intervention_beta_multiplier is not None:
            beta_t = beta_t * cfg.intervention_beta_multiplier(state, t)
        return beta_t * state["S"] * state["I"] / n

    def progression_rate(state: dict, t: float) -> float:
        return cfg.sigma * state["E"]

    def recovery_rate(state: dict, t: float) -> float:
        return cfg.gamma * state["I"]

    reactions = [
        Reaction("infection", infection_rate, {"S": -1, "E": +1}),
        Reaction("progression", progression_rate, {"E": -1, "I": +1}),
        Reaction("recovery", recovery_rate, {"I": -1, "R": +1}),
    ]
    reactions.extend(cfg.extra_reactions)
    return reactions


class StochasticSEIRModel(OutbreakModel):
    """Gillespie-exact stochastic SEIR (or SEIR-derived) model."""

    def __init__(self, config: SEIRConfig, reactions: list[Reaction] | None = None):
        if config.population_size <= 0:
            raise ValueError("population_size must be > 0")
        total_initial = config.initial_exposed + config.initial_infectious
        if total_initial > config.population_size:
            raise ValueError("initial_exposed + initial_infectious cannot exceed population_size")
        self.config = config
        self.reactions = reactions if reactions is not None else build_seir_reactions(config)

    @property
    def population_size(self) -> int:
        return self.config.population_size

    def run(self, rng: np.random.Generator) -> SimulationResult:
        cfg = self.config
        state = {
            "S": cfg.population_size - cfg.initial_exposed - cfg.initial_infectious,
            "E": cfg.initial_exposed,
            "I": cfg.initial_infectious,
            "R": 0,
        }
        # extend state with any compartments introduced by extra_reactions' deltas
        for r in self.reactions:
            for key in r.delta:
                state.setdefault(key, 0)

        t = 0.0
        n_events = 0
        time_series: list[tuple[float, dict]] = [(0.0, dict(state))]
        new_infections_log: list[float] = []  # times at which an "infection" reaction fired

        while t < cfg.max_time and n_events < cfg.max_events:
            rates = np.array([max(r.rate_fn(state, t), 0.0) for r in self.reactions])
            total_rate = rates.sum()
            if total_rate <= 0:
                break  # no reaction can fire -> system at absorbing state

            dt = rng.exponential(1.0 / total_rate)
            t += dt
            if t >= cfg.max_time:
                break

            chosen = rng.choice(len(self.reactions), p=rates / total_rate)
            reaction = self.reactions[chosen]
            for key, d in reaction.delta.items():
                state[key] = state.get(key, 0) + d
            if reaction.name == "infection":
                new_infections_log.append(t)

            n_events += 1
            time_series.append((t, dict(state)))

        final_size = cfg.population_size - state.get("S", 0)
        attack_rate = final_size / cfg.population_size

        if new_infections_log:
            day_bins = np.arange(0, np.ceil(max(new_infections_log)) + 2)
            daily_incidence, _ = np.histogram(new_infections_log, bins=day_bins)
            peak_idx = int(np.argmax(daily_incidence))
            peak_incidence = int(daily_incidence[peak_idx])
            peak_time = float(peak_idx)
        else:
            daily_incidence = np.array([cfg.initial_infectious])
            peak_incidence = cfg.initial_infectious
            peak_time = 0.0

        extinct = state.get("I", 0) == 0 and state.get("E", 0) == 0

        return SimulationResult(
            final_size=final_size,
            attack_rate=attack_rate,
            peak_incidence=peak_incidence,
            peak_time=peak_time,
            duration=t,
            daily_incidence=daily_incidence,
            extinct=extinct,
            generation_sizes=None,
            metadata={"model": "stochastic_seir", "n_events": n_events, "final_state": state},
        )
