"""
Abstract base interface for outbreak models.

Both the stochastic branching process model and the stochastic SEIR model
implement this interface, which lets simulations/monte_carlo.py drive either
one identically. This is the key abstraction that makes the Monte Carlo
engine model-agnostic.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np


@dataclass
class SimulationResult:
    """The output of a single stochastic simulation run (one Monte Carlo iteration).

    Kept intentionally lightweight (arrays, not a full trajectory object with
    behavior) so that thousands of these can be held in memory during a
    Monte Carlo batch without excessive overhead.
    """

    final_size: int  # total number infected over the whole simulation
    attack_rate: float  # final_size / population_size
    peak_incidence: int  # largest number of new infections in a single generation/timestep
    peak_time: float  # generation index or day at which peak_incidence occurred
    duration: float  # generations or days until no infectious individuals remain
    daily_incidence: np.ndarray  # new infections per generation/timestep (calendar-time-mapped for SEIR; generation-indexed for branching)
    extinct: bool  # True if the chain/epidemic died out (as opposed to hitting a size cap)
    generation_sizes: np.ndarray | None = None  # branching-process-specific: cases per generation
    metadata: dict = field(default_factory=dict)


class OutbreakModel(ABC):
    """Common interface for all outbreak transmission models in this project."""

    @abstractmethod
    def run(self, rng: np.random.Generator) -> SimulationResult:
        """Execute a single stochastic realization and return its result.

        Implementations must be fully driven by the passed-in Generator --
        no implementation may call np.random.* directly or construct its own
        Generator internally, since reproducibility (docs/reproducibility.md)
        depends on the caller controlling every source of randomness.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def population_size(self) -> int:
        raise NotImplementedError
