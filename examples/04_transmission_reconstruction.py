"""
Example 4: Reconstruct a probabilistic transmission network from case-timing
data alone (no genomic data required -- see transmission_inference/epi_reconstruction.py),
identify likely superspreaders, and visualize the network.

Run with:  python examples/04_transmission_reconstruction.py
"""

import numpy as np

from outbreak_simulator.data import get_parameter
from outbreak_simulator.transmission_inference import (
    reconstruction_accuracy,
    simulate_ground_truth_tree,
    wallinga_teunis_reconstruction,
)
from outbreak_simulator.visualization import plot_transmission_network

if __name__ == "__main__":
    gen_interval = get_parameter("sars_cov_2", "generation_interval")
    # Gamma method-of-moments shape/scale from the evidence-table mean/CI (same
    # math as models/distributions.py:gamma_from_mean_ci, inlined here since we
    # need the raw (shape, scale) numbers rather than a frozen distribution object)
    import scipy.stats as stats
    z = stats.norm.ppf(0.5 + 0.95 / 2)
    sd = (gen_interval.high - gen_interval.low) / (2 * z)
    shape = gen_interval.point_estimate**2 / sd**2
    scale = sd**2 / gen_interval.point_estimate

    rng = np.random.default_rng(2026)
    case_times, true_infector = None, None
    for seed in range(200):
        rng = np.random.default_rng(seed)
        case_times, true_infector = simulate_ground_truth_tree(
            r_effective=2.5, k_dispersion=0.15,  # SARS-CoV-2-like superspreading dynamics
            generation_interval_shape=shape, generation_interval_scale=scale,
            population_size=80, rng=rng,
        )
        if 15 <= len(case_times) <= 60:
            break

    print(f"Simulated outbreak: {len(case_times)} cases (this is SYNTHETIC data with known ground truth,")
    print("used here to demonstrate and validate the reconstruction method -- see")
    print("docs/validation_plan.md for real accuracy figures measured this way.")
    print()

    reconstructed = wallinga_teunis_reconstruction(case_times, shape, scale)
    acc = reconstruction_accuracy(reconstructed, true_infector)
    print(f"Reconstruction accuracy vs known ground truth: {acc}")

    top_spreaders = sorted(reconstructed.superspreader_scores.items(), key=lambda x: -x[1])[:5]
    print(f"\nTop 5 reconstructed superspreaders (case_id, expected onward transmissions): {top_spreaders}")

    fig = plot_transmission_network(reconstructed, title="Reconstructed transmission network (timing-only)")
    fig.savefig("transmission_network.png", dpi=120)
    print("\nSaved transmission_network.png")

    print("\n--- For real genomic data, see transmission_inference/transphylo_interface.py ---")
    print("(requires R + TransPhylo installed; documented but not exercised by this example)")
