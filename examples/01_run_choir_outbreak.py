"""
Example 1: Run a bundled outbreak scenario with full Monte Carlo uncertainty
propagation, and compare the result against its real-world benchmark.

Run with:  python examples/01_run_choir_outbreak.py
"""

from outbreak_simulator.simulations import print_summary, run_scenario
from outbreak_simulator.validation import calibrate_scenario, print_calibration_report

if __name__ == "__main__":
    result = run_scenario("choir_rehearsal", n_iterations=10_000, seed=20260720)
    print(print_summary(result))
    print()
    print(print_calibration_report(calibrate_scenario(result)))

    print()
    print("--- Trying a different scenario: measles in an under-vaccinated school ---")
    print("Note: this scenario requires you to set susceptible_fraction explicitly")
    print("(see scenarios.yaml assumptions -- population vaccination coverage is a")
    print("policy-relevant input this project deliberately does not hardcode).")
    print()
    measles_result = run_scenario(
        "measles_school_outbreak", n_iterations=10_000, seed=20260720,
        susceptible_fraction=0.08,  # e.g. an 8% under-vaccinated pocket in a school of 400
    )
    print(print_summary(measles_result))
