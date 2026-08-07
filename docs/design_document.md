# Scientific Design Document

## Monte Carlo Infectious Disease Outbreak Simulator

**Status:** v0.1.0 — a working core system covering all required capability
areas, with an intentionally scoped (not exhaustive) literature base. See
"What this project is and isn't" below before using this for anything
consequential.

---

## 1. Purpose

Estimate transmission risk, quantify uncertainty, evaluate interventions,
and (where possible) reconstruct transmission chains for **closed-setting,
introduction-driven outbreaks** — a single or few index cases entering a
defined population (a choir, a barracks, a dormitory, a school), rather than
open, ongoing community transmission. This scope is deliberate: it matches
both the scenarios in the project brief and the settings where a
individual-level, superspreading-aware model earns its complexity relative
to a simpler compartmental model.

## 2. What this project is and isn't

**Is:** A research/educational tool with a real, sourced, citable parameter
base; a working, tested Monte Carlo engine; and an honest accounting of
where the evidence is strong versus thin.

**Isn't:** A validated operational public-health decision-support tool. See
`docs/validation_plan.md` for the specific validation tiers this project
does and does not satisfy. Do not use this to make real intervention
decisions without expert epidemiological review — several parameters
(flagged `low_confidence` throughout `data/parameters/pathogens.yaml`) are
explicitly documented placeholders, not measured values.

## 3. Scientific scope: the six pathogens and eight scenarios

| Pathogen | Scenarios using it | Why this pathogen |
|---|---|---|
| SARS-CoV-2 (ancestral) | Choir rehearsal, military barracks, university dormitory | Best-documented superspreading case studies of any pathogen in this set; three genuinely different settings for the *same* pathogen isolates the effect of setting from the effect of biology. |
| Mpox (2022 clade IIb) | Gathering-associated (sexual network) | A pathogen whose real-world 2022-2026 outbreak was structurally *about* small-gathering/network transmission — the least "generic" scenario in the set. |
| Influenza | Institutional outbreak | Textbook contrast case: low overdispersion (k≈1) vs. SARS-CoV-2's k≈0.15 — same modeling machinery, qualitatively different dynamics. |
| Norovirus | Institutional (school/camp) outbreak | Shortest incubation period, environmentally-transmissible confound (point-source contamination vs. person-to-person) — stress-tests the "single index case" assumption other scenarios rely on. |
| Measles | Under-vaccinated school/community | Highest R0 of any pathogen studied; the vaccination-coverage-as-explicit-input design (see §5) exists largely because of this scenario. |
| Varicella | School/childcare outbreak | High R0, very well-documented household SAR, but a documented closed-setting *counterexample* (military training centre SAR of 21% vs. household 85%) that motivates keeping population-susceptibility separate from setting/contact-structure — see `data/scenarios/scenarios.yaml`. |

For each: required parameters, literature-derived ranges, uncertainty
sources, and assumptions/limitations are all in
`data/parameters/pathogens.yaml` and `data/scenarios/scenarios.yaml` —
**not** duplicated here, because a second copy of the evidence table would
inevitably drift out of sync with the first. This document explains
*rationale*; the YAML files are the source of truth for *values*.

## 4. Literature review methodology (read before trusting any parameter)

This project ran **targeted web searches against primary literature and
public-health-agency guidance**, not a PRISMA-compliant systematic review
(multi-database search, deduplication, two-stage screening). That distinction
matters and is stated in every pathogen's `review_method` field. A genuine
systematic review is a multi-week undertaking for a human team; claiming to
have done one would be dishonest. What *was* done, for every pathogen: at
least one real, retrievable, named source per required parameter, read and
cited, with disagreements between sources recorded rather than papered over
(see e.g. mpox's `infectious_period`, where two sources disagree by a factor
of five, and both are reported).

Every parameter carries an `evidence_quality` tag
(`meta_analysis` / `multi_study` / `single_study` / `expert_consensus` /
`modeled_estimate` / `low_confidence`) precisely so a reader never has to
guess how much weight a number deserves. Nine parameters across the six
pathogens are tagged `low_confidence` — mostly `k_dispersion` for pathogens
without a dedicated overdispersion study (norovirus, measles, varicella,
mpox) and norovirus's community R0. These are the first places to improve
this project with better evidence, not blemishes to hide.

## 5. Key modeling decisions and their rationale

### 5.1 Two-level (outer/inner) Monte Carlo design
Every simulation separates **epistemic uncertainty** ("we don't know the
true R0 precisely") from **aleatory stochasticity** ("even at a known R0,
outcomes vary run to run"). The outer loop draws one value of each uncertain
parameter per Monte Carlo iteration from its evidence-table distribution;
the inner loop runs one stochastic transmission realization at those fixed
values. Conflating these — e.g. reporting simulation-to-simulation variance
as if it were "uncertainty in R0" — is a common and consequential error this
design avoids by construction. See `models/distributions.py` and
`simulations/monte_carlo.py` docstrings for the full argument.

### 5.2 Branching process as the primary model, not SEIR
For the closed, small-population, introduction-driven settings this project
targets, individual-level heterogeneity (superspreading) is often the whole
story — the Skagit choir outbreak's 86% attack rate and the ~35% probability
this project's model assigns to "the chain fizzles at 1 case" both come from
the *same* k=0.15 dispersion parameter. A mean-field SEIR model averages
this away. The stochastic SEIR model (`models/seir.py`) is included and
fully implemented — it is the right tool for larger, longer-duration, or
less superspreading-driven analyses — but the branching process
(`models/branching_process.py`) is what every bundled scenario actually
uses. See `docs/architecture.md` §4 for the full tradeoff discussion,
including why this project does *not* additionally implement a full
agent-based/network model.

### 5.3 Setting (scenario) parameters are separate from pathogen (biological) parameters
`contact_rate_multiplier` on a scenario, not a pathogen field, represents
"how much does *this setting* elevate transmission relative to baseline
community spread." This is what lets `sars_cov_2` be one parameter set
feeding three scenarios (choir/barracks/dormitory) with wildly different
effective R, and is directly motivated by real evidence: norovirus's
person-to-person effective R in a closed school outbreak was measured at
8.92 (Xu et al. 2023) versus a general community R0 in the low single
digits — the *setting* did that, not a different pathogen.

### 5.4 Statistical distribution choices
| Quantity type | Distribution | Why |
|---|---|---|
| Bounded probability (SAR, efficacy, coverage) | Beta | Natural [0,1] support; two parameters fit a mean + CI width without extra degrees of freedom. |
| Duration: incubation period | Lognormal | Matches the convention used directly by the cited primary literature (e.g. Lauer et al.-style COVID incubation studies). |
| Duration: infectious/generation/serial interval | Gamma | Matches convention in the cited generation-interval literature (e.g. Ganyani et al.). |
| Individual-level offspring (secondary case) count | Negative Binomial | The standard mechanism for superspreading (Lloyd-Smith et al. 2005): mean=R, dispersion=k, k→∞ recovers Poisson. |
| Compartmental transition counts (SEIR) | Poisson (via Gillespie SSA) | Homogeneous-mixing limit — appropriate once individual identity no longer matters, which is what distinguishes the SEIR pathway from the branching process pathway. |

### 5.5 Reproducibility architecture
Every Monte Carlo iteration's randomness derives from a single master seed
via `numpy.random.SeedSequence.spawn()` — the modern numpy-recommended
pattern, superseding global `np.random.seed()`. Two independent child
streams are spawned per iteration (one for outer-loop parameter sampling,
one for inner-loop stochastic simulation) so that reproducibility doesn't
depend on call-order coincidences. See `docs/reproducibility.md`.

## 6. Explicit non-goals

- **Not** a real-time forecasting tool (no calendar-time-anchored nowcasting).
- **Not** a spatial/geographic model (no travel, no cross-setting mixing).
- **Not** a within-host/immunological model.
- **Not** validated for operational public-health decision-making (see
  `docs/validation_plan.md`).
- **Does not** claim the mpox `k_dispersion`, or the measles/varicella/
  norovirus `k_dispersion` values, are literature-derived — they are
  explicitly flagged, reasoned placeholders (see §4).
