"""Validation: goodness-of-fit metrics and tiered scenario calibration against real-world observed outcomes."""

from outbreak_simulator.validation.calibration import (
    ObservationCalibration,
    ScenarioCalibrationReport,
    calibrate_scenario,
    print_calibration_report,
)
from outbreak_simulator.validation.metrics import (
    CoverageResult,
    posterior_predictive_check,
    predictive_coverage,
    relative_bias,
    root_mean_squared_error,
)

__all__ = [
    "ObservationCalibration", "ScenarioCalibrationReport", "calibrate_scenario", "print_calibration_report",
    "CoverageResult", "posterior_predictive_check", "predictive_coverage", "relative_bias", "root_mean_squared_error",
]
