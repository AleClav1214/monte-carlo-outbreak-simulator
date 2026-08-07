"""
Epidemiological transmission-chain reconstruction (no genomic data required).

Implements the Wallinga & Teunis (2004, Am J Epidemiol, "Different
epidemic curves for severe acute respiratory syndrome reveal similar
impacts of control measures") method: given only case symptom-onset times
and a generation/serial-interval distribution, compute the relative
likelihood that case j was infected by case i as proportional to the
serial-interval density evaluated at (t_j - t_i), normalized across all
candidate infectors of j. This produces a probabilistic "who-infected-whom"
network without any genomic sequence data -- the required
requirement-#10 alternative for when genomic data is unavailable.

This is directly testable against this project's own branching-process
simulator (which records each case's true infector and infection time),
letting us check reconstruction accuracy against known ground truth --
something that is NOT possible for the TransPhylo/genomic pathway, since
this project has no real genomic outbreak dataset (see
transphylo_interface.py for that documented-but-unvalidated interface).
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
from scipy import stats


@dataclass
class ReconstructedTransmissionNetwork:
    graph: nx.DiGraph  # nodes = case indices, edge (i,j) weight = P(i infected j)
    most_likely_infector: dict[int, int | None]  # case_id -> most probable infector (None for presumed index case(s))
    superspreader_scores: dict[int, float]  # case_id -> expected number of onward transmissions attributed to them


def wallinga_teunis_reconstruction(
    case_times: np.ndarray,
    generation_interval_shape: float,
    generation_interval_scale: float,
) -> ReconstructedTransmissionNetwork:
    """Reconstruct a probabilistic transmission network from case timing data alone.

    Parameters
    ----------
    case_times: array of symptom-onset (or infection) times, one per case,
        in the same time unit as the generation-interval parameters.
    generation_interval_shape, generation_interval_scale: Gamma distribution
        parameters for the serial/generation interval (this project's
        pathogen evidence tables report generation_interval or
        serial_interval as Gamma-distributed -- see data/parameters/pathogens.yaml).

    Method
    ------
    For every ordered pair (i, j) with t_i < t_j, compute
        w_ij = Gamma_pdf(t_j - t_i; shape, scale)
    Then normalize across all candidates for each j:
        p_ij = w_ij / sum_i' w_i'j
    This is exactly the Wallinga-Teunis relative-likelihood construction.
    Cases with no earlier case within a plausible generation-interval window
    are treated as (co-)index cases with no assigned infector.
    """
    n = len(case_times)
    order = np.argsort(case_times)
    sorted_times = case_times[order]

    gen_dist = stats.gamma(a=generation_interval_shape, scale=generation_interval_scale)

    graph = nx.DiGraph()
    for idx in range(n):
        graph.add_node(int(order[idx]), time=float(sorted_times[idx]))

    most_likely_infector: dict[int, int | None] = {}
    for j in range(n):
        t_j = sorted_times[j]
        candidates = [i for i in range(j) if sorted_times[i] < t_j]
        if not candidates:
            most_likely_infector[int(order[j])] = None
            continue
        weights = np.array([gen_dist.pdf(t_j - sorted_times[i]) for i in candidates])
        total = weights.sum()
        if total <= 0:
            most_likely_infector[int(order[j])] = None
            continue
        probs = weights / total
        for c_idx, p in zip(candidates, probs):
            if p > 1e-6:
                graph.add_edge(int(order[c_idx]), int(order[j]), weight=float(p))
        best = candidates[int(np.argmax(probs))]
        most_likely_infector[int(order[j])] = int(order[best])

    superspreader_scores: dict[int, float] = {int(node): 0.0 for node in graph.nodes}
    for u, v, data in graph.edges(data=True):
        superspreader_scores[u] = superspreader_scores.get(u, 0.0) + data["weight"]

    return ReconstructedTransmissionNetwork(
        graph=graph, most_likely_infector=most_likely_infector, superspreader_scores=superspreader_scores,
    )


def simulate_ground_truth_tree(
    r_effective: float,
    k_dispersion: float,
    generation_interval_shape: float,
    generation_interval_scale: float,
    population_size: int,
    rng: np.random.Generator,
    max_generations: int = 50,
) -> tuple[np.ndarray, dict[int, int | None]]:
    """Individually-tracked branching-process simulation used ONLY to
    generate known-ground-truth (case_times, true_infector) pairs for
    testing wallinga_teunis_reconstruction()'s accuracy (see
    reconstruction_accuracy() below and tests/unit/test_transmission_inference.py).

    This intentionally duplicates the transmission mechanics of
    models/branching_process.py (same NB-offspring + susceptible-depletion
    logic) but additionally records each individual case's specific
    infector, which the main model does not track (by design -- it operates
    on generation-aggregate counts for performance, since individual
    identity is not needed for the Monte Carlo attack-rate/uncertainty
    outputs that are that model's actual purpose). Keeping this tracked
    variant separate avoids adding per-individual bookkeeping overhead to
    the hot path of the main Monte Carlo engine.
    """
    from outbreak_simulator.models.distributions import negative_binomial_offspring

    case_times = [0.0]
    true_infector: dict[int, int | None] = {0: None}
    frontier = [0]  # indices of currently-infectious individuals
    susceptible = population_size - 1
    next_id = 1

    for _ in range(max_generations):
        if not frontier or susceptible <= 0:
            break
        depletion = susceptible / population_size
        r_eff_gen = r_effective * depletion
        offspring_counts = negative_binomial_offspring(r_eff_gen, k_dispersion, rng, size=len(frontier))

        new_frontier = []
        for infector_idx, n_offspring in zip(frontier, offspring_counts):
            n_actual = min(int(n_offspring), susceptible)
            for _ in range(n_actual):
                delay = rng.gamma(generation_interval_shape, generation_interval_scale)
                new_time = case_times[infector_idx] + delay
                case_times.append(float(new_time))
                true_infector[next_id] = infector_idx
                new_frontier.append(next_id)
                next_id += 1
                susceptible -= 1
                if susceptible <= 0:
                    break
            if susceptible <= 0:
                break
        frontier = new_frontier

    return np.array(case_times), true_infector


def reconstruction_accuracy(
    reconstructed: ReconstructedTransmissionNetwork, true_infector: dict[int, int | None]
) -> dict:
    """Compare reconstructed most-likely infectors against KNOWN ground truth
    (available only from simulated data, e.g. models/branching_process.py
    output, where the true transmission tree is recorded by construction --
    this is how this project validates the reconstruction method itself,
    distinct from validating outbreak-size predictions).
    """
    matches = 0
    total = 0
    for case_id, true_inf in true_infector.items():
        if true_inf is None:
            continue  # skip true index cases; reconstruction correctly has no infector to guess for them by construction
        total += 1
        if reconstructed.most_likely_infector.get(case_id) == true_inf:
            matches += 1
    accuracy = matches / total if total > 0 else float("nan")
    return {"n_evaluated": total, "n_correct": matches, "accuracy": accuracy}
