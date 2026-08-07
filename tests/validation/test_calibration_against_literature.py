"""
Validation tests.

These tests check the project's outputs against real, literature-reported
benchmarks (scenarios.yaml observed_outcomes) and against known ground
truth (for the transmission-reconstruction method, generated synthetically
since no real genomic/timing dataset with confirmed ground truth is
available -- see transmission_inference/epi_reconstruction.py docstring).

Per this project's validation-tier discipline: passing these tests
demonstrates the model is not GROSSLY miscalibrated against the available
real-world benchmarks. It does NOT demonstrate general external validity
(see docs/validation_plan.md) -- these are consistency checks, not proof of
correctness.
"""

from __future__ import annotations

import numpy as np

from outbreak_simulator.simulations import run_scenario
from outbreak_simulator.transmission_inference import (
    reconstruction_accuracy,
    simulate_ground_truth_tree,
    wallinga_teunis_reconstruction,
)
from outbreak_simulator.validation import calibrate_scenario


class TestCalibrationAgainstLiteratureBenchmarks:
    """For each scenario with a real-world benchmark, check the observed
    value is not a wild outlier relative to the model's predictive
    distribution (i.e. PPC p-value not vanishingly small) -- a coarse but
    meaningful 'not grossly wrong' check."""

    def _check_scenario(self, scenario_id: str, min_ppc_pvalue: float = 0.01):
        result = run_scenario(scenario_id, n_iterations=5000, seed=2026)
        report = calibrate_scenario(result)
        assert report.n_independent_benchmarks >= 1
        for obs in report.per_observation:
            assert obs.ppc_pvalue >= min_ppc_pvalue, (
                f"{scenario_id}: observed value {obs.coverage.observed_value:.1%} has PPC p={obs.ppc_pvalue:.4f}, "
                f"suggesting it is a poor fit to the model's predictive distribution -- investigate before trusting "
                f"this scenario's parameters."
            )

    def test_choir_rehearsal_calibration(self):
        self._check_scenario("choir_rehearsal")

    def test_military_barracks_calibration(self):
        self._check_scenario("military_barracks")

    def test_university_dormitory_calibration(self):
        self._check_scenario("university_dormitory")

    def test_norovirus_outbreak_calibration(self):
        self._check_scenario("norovirus_outbreak")

    def test_varicella_school_outbreak_calibration(self):
        self._check_scenario("varicella_school_outbreak")


class TestTransmissionReconstructionAccuracy:
    """Checks the epidemiological (non-genomic) reconstruction method
    against known ground truth from synthetic outbreaks. Documents the
    REAL, measured performance rather than an assumed one -- see module
    docstring in transmission_inference/epi_reconstruction.py for why
    accuracy is expected to be well above chance but well below perfect."""

    def test_reconstruction_beats_random_baseline_at_realistic_scale(self):
        accuracies = []
        for seed in range(20):
            rng = np.random.default_rng(seed + 5000)
            case_times, true_infector = simulate_ground_truth_tree(
                r_effective=2.0, k_dispersion=0.5,
                generation_interval_shape=4.0, generation_interval_scale=1.3,
                population_size=60, rng=rng,
            )
            if len(case_times) < 15:
                continue
            reconstructed = wallinga_teunis_reconstruction(case_times, 4.0, 1.3)
            acc = reconstruction_accuracy(reconstructed, true_infector)
            if acc["n_evaluated"] > 0:
                accuracies.append(acc["accuracy"])

        assert len(accuracies) >= 3, "not enough non-trivial outbreaks generated across seeds to evaluate"
        mean_accuracy = float(np.mean(accuracies))
        # Random-guess baseline is well below this at realistic candidate-set sizes;
        # this asserts *meaningfully above chance*, not near-perfect reconstruction,
        # which matches the real literature on timing-only (non-genomic) reconstruction.
        assert mean_accuracy > 0.15, (
            f"mean reconstruction accuracy {mean_accuracy:.2f} is at or below what a naive random guess "
            f"among plausible candidates would achieve -- the method may be broken, not just imprecise"
        )

    def test_index_case_correctly_has_no_assigned_infector(self):
        rng = np.random.default_rng(1)
        case_times, true_infector = simulate_ground_truth_tree(2.0, 0.5, 4.0, 1.3, 60, rng)
        reconstructed = wallinga_teunis_reconstruction(case_times, 4.0, 1.3)
        assert reconstructed.most_likely_infector[0] is None
