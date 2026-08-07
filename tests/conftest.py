"""Shared pytest fixtures."""

from __future__ import annotations

import numpy as np
import pytest

from outbreak_simulator.data import get_pathogen, get_scenario


@pytest.fixture
def rng():
    return np.random.default_rng(20260720)


@pytest.fixture
def sars_cov_2():
    return get_pathogen("sars_cov_2")


@pytest.fixture
def measles():
    return get_pathogen("measles")


@pytest.fixture
def choir_scenario():
    return get_scenario("choir_rehearsal")


@pytest.fixture
def barracks_scenario():
    return get_scenario("military_barracks")
