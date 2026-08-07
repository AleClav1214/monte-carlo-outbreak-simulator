"""
Stochastic branching process model (individual-level transmission).

This is the primary model for the superspreading-focused scenarios in this
project (choir rehearsal, mpox gathering, etc.), because it represents
transmission heterogeneity directly via the negative-binomial offspring
distribution (Lloyd-Smith et al. 2005) rather than averaging it away, as a
mean-field compartmental model would.

Algorithm (generation-based branching process with susceptible depletion)
---------------------------------------------------------------------------
1. Generation 0 = the scenario's initial_cases.
2. For each infectious individual in the current generation, draw a
   *potential* offspring count from NegativeBinomial(mean=R_effective,
   dispersion=k). R_effective = R_scenario * (S / N), i.e. scaled down as
   the pool of remaining susceptibles shrinks -- the same mean-field
   depletion correction used in compartmental models, applied here at the
   individual level. This is a standard simplification (see e.g. Kucharski
   et al. 2020 Lancet Infect Dis, who use an equivalent construction) and
   is documented as such rather than presented as a first-principles
   contact-network model.
3. Sum potential offspring across all infectors in the generation; cap the
   *actual* new-case count at the remaining susceptible population (extra
   potential transmissions "hit" already-infected/immune individuals and
   produce no new case) -- the discrete-time analogue of the Reed-Frost
   chain-binomial mechanism.
4. Advance one generation; repeat until either no new cases occur
   (extinction) or the susceptible pool is exhausted or max_generations is
   reached (safety cap).
5. Each case is optionally assigned a calendar time by drawing a
   generation-interval delay from its infector's time -- this produces the
   daily_incidence curve used for epidemic-curve visualization, converting
   the model's natural generation-indexed output into calendar time.

Interventions are applied as multiplicative modifiers to R_effective and/or
k at simulation-construction time (see interventions/ for how modifiers are
composed) -- the branching process model itself has no knowledge of what an
"intervention" is; it just consumes an already-composed effective R and k.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from outbreak_simulator.models.base import OutbreakModel, SimulationResult
from outbreak_simulator.models.distributions import negative_binomial_offspring


@dataclass
class BranchingProcessConfig:
    population_size: int
    initial_cases: int
    r_effective: float  # already composed: baseline_r0 * contact_rate_multiplier * intervention_modifiers
    k_dispersion: float
    generation_interval_sampler: callable | None = None  # rng -> float days; if None, generation index is used as "time"
    max_generations: int = 200
    deplete_susceptibles: bool = True


class BranchingProcessModel(OutbreakModel):
    """Stochastic, individual-level, generation-based branching process with
    negative-binomial offspring (superspreading) and susceptible depletion."""

    def __init__(self, config: BranchingProcessConfig):
        if config.population_size <= 0:
            raise ValueError("population_size must be > 0")
        if config.initial_cases <= 0 or config.initial_cases > config.population_size:
            raise ValueError("initial_cases must be in [1, population_size]")
        if config.k_dispersion <= 0:
            raise ValueError("k_dispersion must be > 0")
        if config.r_effective < 0:
            raise ValueError("r_effective must be >= 0")
        self.config = config

    @property
    def population_size(self) -> int:
        return self.config.population_size

    def run(self, rng: np.random.Generator) -> SimulationResult:
        cfg = self.config
        n = cfg.population_size
        susceptible = n - cfg.initial_cases
        current_generation_size = cfg.initial_cases
        generation_sizes = [current_generation_size]
        case_times = [0.0] * cfg.initial_cases  # calendar time of each case, generation 0 = t=0
        frontier_times = list(case_times)  # times of the individuals currently able to transmit

        total_cases = cfg.initial_cases
        peak_incidence = current_generation_size
        peak_time = 0.0
        extinct = False

        for gen in range(1, cfg.max_generations + 1):
            if current_generation_size == 0 or susceptible <= 0:
                extinct = susceptible > 0 or current_generation_size == 0
                break

            depletion_factor = (susceptible / n) if cfg.deplete_susceptibles else 1.0
            r_eff_this_gen = cfg.r_effective * depletion_factor

            potential_offspring = negative_binomial_offspring(
                r_eff_this_gen, cfg.k_dispersion, rng, size=current_generation_size
            )
            total_potential = int(potential_offspring.sum())
            actual_new_cases = min(total_potential, susceptible)

            if actual_new_cases <= 0:
                extinct = True
                break

            # assign calendar times to the new cases: each new case's time =
            # (a randomly chosen infector's time from the previous frontier) +
            # a generation-interval delay
            new_times: list[float] = []
            if cfg.generation_interval_sampler is not None and frontier_times:
                infector_indices = rng.integers(0, len(frontier_times), size=actual_new_cases)
                delays = np.array([cfg.generation_interval_sampler(rng) for _ in range(actual_new_cases)])
                new_times = list(np.array(frontier_times)[infector_indices] + delays)
            else:
                new_times = [float(gen)] * actual_new_cases  # fall back to generation index as "time"

            susceptible -= actual_new_cases
            total_cases += actual_new_cases
            generation_sizes.append(actual_new_cases)
            frontier_times = new_times
            case_times.extend(new_times)

            if actual_new_cases > peak_incidence:
                peak_incidence = actual_new_cases
                peak_time = float(np.mean(new_times)) if new_times else float(gen)

            current_generation_size = actual_new_cases
        else:
            # loop completed without break -> hit max_generations without dying out
            extinct = False

        case_times_arr = np.array(case_times)
        if cfg.generation_interval_sampler is not None and len(case_times_arr) > 0:
            max_day = int(np.ceil(case_times_arr.max())) + 1
            daily_incidence = np.bincount(np.clip(case_times_arr, 0, None).astype(int), minlength=max_day)
        else:
            daily_incidence = np.array(generation_sizes)

        return SimulationResult(
            final_size=total_cases,
            attack_rate=total_cases / n,
            peak_incidence=peak_incidence,
            peak_time=peak_time,
            duration=float(len(generation_sizes) - 1),
            daily_incidence=daily_incidence,
            extinct=extinct,
            generation_sizes=np.array(generation_sizes),
            metadata={"model": "branching_process", "n_generations": len(generation_sizes) - 1},
        )
