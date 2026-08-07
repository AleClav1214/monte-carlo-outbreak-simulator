# Validation Plan

This document follows a **tiered validation framework**: every validation
claim is classified as `necessary`, `recommended`, `optional`, or `not
currently justified` for the claims this project actually makes, rather
than asserting a blanket "validated" label. A model can be legitimately
useful and honestly documented without being externally validated in the
strong sense — the goal here is precision about which is true where.

## 1. Validation tiers applied to this project

| Tier | Status | What it checks | Where |
|---|---|---|---|
| **Internal resampling** | ✅ Done | Do the model's own mechanics behave as designed (NB offspring mean/variance match theory, population conservation, reproducibility under fixed seeds)? | `tests/unit/test_models.py`, `tests/unit/test_monte_carlo.py` |
| **Internal holdout / recovery** | ✅ Done | Given known parameters, does the model recover the expected qualitative behavior (superspreading signature, extinction probability, saturation at high R)? | `tests/unit/test_models.py::TestBranchingProcessModel` |
| **External cohort validation** | ⚠️ Partial, n=1-2 per scenario | Does the model's predictive distribution contain the real, literature-reported outcome for this exact setting? | `validation/calibration.py`, `tests/validation/test_calibration_against_literature.py` |
| **Site-based validation** | ❌ Not currently justified | Does calibration hold across *multiple independent occurrences* of the same scenario type? | Most scenarios have only one documented real-world instance (see §3). |
| **Orthogonal validation** | ❌ Not attempted | Does this model agree with other established outbreak-modeling tools/papers on the same scenario? | Not in scope for this version. |
| **Temporal validation** | ❌ Not attempted | Does calibration hold across pathogen variants / time periods / population-immunity states? | All SARS-CoV-2 parameters are ancestral-strain 2020 defaults; not tested against later variants. |
| **Functional validation** | ❌ Not currently justified | Are outputs reliable enough to inform real decisions? | Out of scope — this is a research/educational tool. |
| **Translational validation** | ❌ Not attempted | Has this been compared against, or reviewed alongside, tools actually used operationally by public health agencies? | Not attempted. |

**Reporting rule applied throughout this project:** `validation/calibration.py`
never emits an unqualified "validated" — every `ScenarioCalibrationReport`
explicitly lists what it does *not* establish (see its `NOT_established`
field), and every report generated states its validation tier and
benchmark count explicitly.

## 2. What "external cohort validation, n=1" actually checked, per scenario

Method: for each scenario, run 5,000-10,000 Monte Carlo iterations, then
check whether the real-world observed attack rate falls within the
simulated distribution's 95% predictive interval, its percentile rank
within that distribution, and a two-sided posterior-predictive-check
p-value (`validation/metrics.py:predictive_coverage`,
`posterior_predictive_check`).

**Real results from running this** (seed=2026, n=5000 per scenario; see
`tests/validation/test_calibration_against_literature.py` for the automated
version of this check, gated at PPC p ≥ 0.01):

- **Choir rehearsal:** observed 86% attack rate falls at approximately the
  69th percentile of the simulated distribution — comfortably central, not
  a tail event.
- **Military barracks:** both observed benchmarks (22.3% Fort Benning,
  12.8% South Korea) land near the 54th percentile — also comfortably
  central.

**An important honest finding, not smoothed over:** the choir scenario's
simulated attack-rate distribution is strongly **bimodal**, not unimodal —
median outcome is near 1/61 (the chain fizzles immediately) while roughly a
quarter of realizations reach full population saturation (attack rate =
1.0), with comparatively little simulated mass at intermediate values like
the observed 86%. This is the correct qualitative signature of low-*k*
superspreading dynamics in a small population (see Althouse et al. 2020 on
stochastic superspreading dynamics for the general phenomenon), but it also
means **"the observed value is within the 95% interval" is a much weaker
statement here than it would be for a unimodal, well-behaved output** — it
doesn't mean the observed value is near the *typical* simulated outcome,
because for this scenario there isn't a single typical outcome. Report
percentile rank and the full distributional shape alongside coverage, never
coverage alone, for scenarios like this one — `print_calibration_report()`
does this by default.

## 3. What is explicitly NOT established

- **n=1 (or n=2) real-world benchmarks cannot distinguish "this setting is
  reliably modeled" from "this specific documented case happened to be
  consistent with these parameters."** This is the single most important
  caveat in this document. A model that is wrong in a way that happens to
  still be consistent with one real data point will pass every check in §2.
- **The choir scenario's `contact_rate_multiplier` was chosen with the
  Skagit outcome in view.** It was not derived from independent
  first-principles aerosol physics and then checked against the outcome —
  it is closer to a calibration target for that one scenario than an
  out-of-sample prediction. This is stated explicitly in
  `scenarios.yaml: choir_rehearsal.limitations` and repeated here because
  it is easy to lose track of when reading only the calibration report.
- **The mpox gathering scenario has no single-event benchmark at all** —
  only a general rate statement (WHO's ~75% sexual-contact SAR pattern
  across many events, not one documented gathering's investigation). Its
  calibration tier is explicitly weaker than the other scenarios' — see
  `scenarios.yaml: mpox_gathering.limitations`.

## 4. Transmission reconstruction accuracy (a genuinely measured result)

`transmission_inference/epi_reconstruction.py`'s Wallinga-Teunis-based
method was validated against **known ground truth from synthetic outbreaks**
(real transmission trees are not published for any of this project's
scenarios, so this is the only ground truth available — see
`simulate_ground_truth_tree()`). Measured results:

- **At large-epidemic scale (~300 cases):** exact most-likely-infector
  match accuracy ≈ **3.7%** — many temporally-overlapping candidate
  infectors make pairwise disambiguation from timing alone genuinely hard,
  consistent with the transmission-inference literature's standard
  motivation for incorporating genomic data.
- **At this project's actual scale (15-55 cases, matching the bundled
  scenarios):** exact match accuracy ≈ **25-38%** across repeated trials —
  meaningfully above a random-guess baseline, but far from perfect.
  `tests/validation/test_calibration_against_literature.py::TestTransmissionReconstructionAccuracy`
  encodes this as a regression test (asserts >15%, i.e. "clearly better
  than chance," not a much stronger and less defensible claim).

This is reported as a measured number, not an assumed capability, precisely
because the module's own docstring warns that timing-only reconstruction is
a known-weaker method than genomic reconstruction — the measurement
confirms that warning rather than contradicting it.

## 5. Goodness-of-fit / calibration metrics implemented

See `validation/metrics.py`: predictive coverage + percentile rank, RMSE
(for multi-scenario point-estimate comparison), relative bias, and a
posterior-predictive-check-style two-sided p-value. Each function's
docstring states what it can and cannot support — e.g. RMSE is documented as
requiring multiple paired scenario/observation points to be meaningful, and
is not used for any single-scenario n=1 comparison in this project.

## 6. Known model limitations (consolidated from per-file notes)

- Susceptible depletion in the branching process model uses a simple
  mean-field `S/N` scaling of R at each generation — a standard
  simplification, not a full contact-network depletion model.
- The Wells-Riley ventilation approximation
  (`interventions/environmental.py`) uses the small-exponent linear
  approximation (`risk ∝ 1/ACH`), which weakens for long/high-risk exposures
  — flagged explicitly in that module's docstring, and most relevant to the
  choir scenario's 2.5-hour exposure.
- Intervention effects combine multiplicatively assuming independence
  (`interventions/stack.py`) — real interventions can interact
  non-independently (e.g. masking matters less once ventilation already
  dilutes aerosol to low concentration).
- No age-structured contact patterns anywhere in the project — relevant
  especially to the measles/varicella school scenarios.
