# Evidence Tables

**Generated from `data/parameters/pathogens.yaml` — do not hand-edit this file.**
Re-run `python scripts/generate_evidence_tables.py` after changing the YAML.

See `data/parameters/pathogens.yaml`'s file header for the full literature-review
methodology statement (targeted search, not a PRISMA systematic review).

## Seasonal influenza A/B

*respiratory virus (droplet + airborne). Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Incubation period | 2 | [1, 4] | days | Expert consensus | Standard public-health guidance (CDC, WHO); well-established, low between-study varianc... |
| Infectious period | 5 | [3, 7] | days | Expert consensus | Standard public-health guidance: infectious from ~1 day before symptom onset to ~5-7 da... |
| Negative-binomial dispersion parameter (k) | 1 | [0.6, 2] | dimensionless | Single study | Comparative overdispersion estimate for 1918 pandemic influenza, k~1 (approximately hom... |
| Basic reproduction number (R0) | 1.3 | [0.9, 2.1] | dimensionless | Meta-analysis | Biggerstaff et al. 2014 BMC Infect Dis systematic review of R estimates for seasonal, p... |
| Household secondary attack rate | 0.1 | [0.05, 0.19] | probability | Multi-study | Multiple household-transmission studies: San Antonio H1N1 (4-13% by case definition), N... |
| Serial interval | 3 | [2.6, 4] | days | Multi-study | Household studies: Ontario H1N1 (median 3.0 days), Navarra H1N1pdm09 (mean 3.7 days, 95... |

**Justifications and notes:**

- **Incubation period:** Textbook-stable figure; lognormal for consistency with other respiratory-virus incubation periods in this project.
- **Infectious period:** Textbook-stable figure.
- **Negative-binomial dispersion parameter (k):** Influenza is the standard textbook contrast case for "low overdispersion" transmission (relative to SARS-like coronaviruses), attributed to influenza's efficient, less contact-structure-dependent respiratory spread. Value used here reflects that qualitative consensus more than a precise, dedicated influenza-k estimation study.
  *Notes:* Secondary source relaying the Lloyd-Smith et al. 2005 comparative figure, not the primary paper directly -- flagged accordingly.
- **Basic reproduction number (R0):** Systematic review is the most authoritative available source for influenza R estimates across many outbreak settings and years.
  *Notes:* Pandemic strains are higher: 2009 H1N1 commonly estimated at 1.4-1.6 (some early estimates up to 2.9); the 1918 pandemic strain had a median R of 1.80 (IQR 1.47-2.27) per the same systematic review. Use a pandemic-strain override for pandemic scenarios.
- **Household secondary attack rate:** Beta distribution for a bounded probability; point estimate and interval triangulated across four independent household-contact studies plus the WHO seasonal benchmark range, deliberately centered on the seasonal (non-pandemic) figure.
  *Notes:* Pandemic-strain household SAR runs meaningfully higher (up to ~33% per WHO H1N1 estimates) -- override for pandemic scenarios.
- **Serial interval:** Gamma distribution, standard for serial-interval/generation-interval quantities; range spans the two cited study point estimates.

## Measles (rubeola)

*respiratory virus (airborne, highly efficient). Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Incubation period (exposure to prodromal symptom onset) | 10 | [8, 12] | days | Expert consensus | Louisiana Dept. of Health epidemiology manual / CDC Pink Book: incubation 8-12 days to ... |
| Infectious period (relative to rash onset) | 8 | [6, 9] | days | Expert consensus | Standard public-health guidance: infectious ~4 days before through ~4 days after rash o... |
| Negative-binomial dispersion parameter (k) | 1 | [0.3, 2] | dimensionless | ⚠️ LOW CONFIDENCE | No dedicated measles k-estimation study located in this search pass |
| Basic reproduction number (R0) | 15 | [12, 18] | dimensionless | Expert consensus | Consistently cited across independent sources: CDC EID (Voigt et al. 2024, PMC11346981)... |
| Household secondary attack rate, susceptible (unvaccinated) contacts | 0.85 | [0.75, 0.9] | probability | Expert consensus | Cross-pathogen SAR benchmarking (measles ~75% commonly cited as a conservative figure; ... |
| Single-dose MMR vaccine efficacy (historical outbreak-derived estimate) | 0.973 | [0.801, 0.999] | probability | Single study | Cherry et al., secondary attack rates during the 1974 Cheyenne/Standing Rock Sioux Rese... |

**Justifications and notes:**

- **Incubation period (exposure to prodromal symptom onset):** Stable, widely-replicated public-health-agency figure.
  *Notes:* IMPORTANT MODELING DETAIL: measles cases are infectious from ~4 days BEFORE rash onset, i.e. before the classic diagnostic sign appears -- this project's infectious_period window is centered on this pre-rash infectious period, not on rash-to-recovery, since that is when most transmission actually happens.
- **Infectious period (relative to rash onset):** Textbook-stable figure (CDC Pink Book, WHO).
- **Negative-binomial dispersion parameter (k):** PLACEHOLDER ASSUMPTION based on general epidemiological reasoning rather than a literature point estimate: pathogens with very high, airborne-efficient baseline R0 that can infect most susceptibles in a room tend to show somewhat less *relative* overdispersion than close-contact-dependent pathogens, because the transmission opportunity is less bottlenecked by individual behavioral heterogeneity. This is an analyst judgment call, not a citation -- flagged accordingly. Documented superspreading events do occur (e.g. single cases in clinics/schools infecting dozens of susceptible contacts), so this should not be read as "measles has no superspreading."
  *Notes:* Highest-priority parameter to replace with a real literature estimate if this project is extended.
- **Basic reproduction number (R0):** One of the most robust, widely-replicated figures in infectious disease epidemiology; multiple independent sources across different decades converge tightly on 12-18.
  *Notes:* Among the highest R0 values of any known human pathogen; the basis for the ~95% vaccination coverage threshold commonly cited for measles herd immunity.
- **Household secondary attack rate, susceptible (unvaccinated) contacts:** Beta distribution for bounded probability; among the most consistently-cited SAR figures across sources given how extensively measles household transmission has been studied historically.
- **Single-dose MMR vaccine efficacy (historical outbreak-derived estimate):** Directly usable as the default vaccination-intervention efficacy parameter for measles scenarios; single historical study, but consistent with the broadly-cited modern two-dose MMR efficacy figure (~97%) from routine immunization program data.
  *Notes:* Two-dose modern MMR efficacy is generally cited even higher (~97-99%); this single-dose 1970s figure is retained as a conservative default.

## Mpox (2022-2024 outbreak clade IIb, human-to-human transmission)

*orthopoxvirus (close/skin contact, respiratory droplet, fomite). Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Incubation period | 8 | [6, 13] | days | Multi-study | CDC Center for Forecasting and Outbreak Analytics household-transmission model (mean 8 ... |
| Infectious period (symptomatic through lesion crusting) | 18 | [14, 27] | days | ⚠️ LOW CONFIDENCE | CDC household-transmission model (mean 27 days); WHO/CDC clinical guidance (infectious ... |
| Negative-binomial dispersion parameter (k) | 0.5 | [0.1, 1.5] | dimensionless | ⚠️ LOW CONFIDENCE | No consensus estimate located in literature search (see Alshahrani et al. 2025 PMC12056... |
| Basic reproduction number (R0) | 2.2 | [1.37, 3.68] | dimensionless | Multi-study | Multi-country 2022 estimates: UK 2.32 (deep-learning approach), Berlin 2.13 (95% PrI 1.... |
| Household (non-sexual close-contact) secondary attack rate | 0.08 | [0, 0.2] | probability | Multi-study | Systematic review of endemic-country household data (pooled 8%, range 0-11%, Guagliardo... |
| Secondary attack rate among sexual contacts (gathering-associated context) | 0.75 | [0.6, 0.85] | probability | Single study | WHO Multi-country External Situation Report (secondary attack rate among sexual contact... |

**Justifications and notes:**

- **Incubation period:** Lognormal for consistency with the general convention for incubation periods (right-skewed). Point estimate and interval bridge the CDC modeling assumption (8 days) and the somewhat higher deep-learning estimate (9-10 days), with the WHO-cited upper tail (13 days for 95% of cases; up to 21 days maximum) as the high bound.
- **Infectious period (symptomatic through lesion crusting):** Genuinely ambiguous parameter: some sources describe a short acute febrile-prodrome infectious window (~4-5 days) while others (and formal isolation guidance) treat the infectious period as lasting until full lesion crusting, i.e. 2-4 weeks. This project uses the longer, precautionary-isolation-guidance-consistent value because it is the more consequential one for outbreak-size projections; a shorter value would be optimistic and is flagged here as an important sensitivity-analysis target.
  *Notes:* FLAGGED FOR SENSITIVITY ANALYSIS: outbreak size projections are likely sensitive to this parameter given the wide literature range. Run sensitivity/one_way.py on this parameter before drawing conclusions from any mpox scenario.
- **Negative-binomial dispersion parameter (k):** PLACEHOLDER ASSUMPTION, explicitly flagged. Set to a moderate-overdispersion value (higher/less extreme than SARS-CoV-2's k~0.1) reasoning by analogy from the sexual-network-concentrated transmission pattern (which typically implies some heterogeneity, since a minority of highly-connected individuals in the network drive disproportionate spread) without a directly measured k. DO NOT treat this as a literature-derived estimate; update as soon as a primary k-estimation study for mpox becomes available.
  *Notes:* This is the single weakest-evidence parameter in the mpox parameter set.
- **Basic reproduction number (R0):** No single pooled meta-analytic R0 exists for the 2022-2024 outbreak; point estimate is the rough center of mass across ~6 independent country/method estimates found. Wide interval deliberately preserves the real between-country heterogeneity (driven by differing sexual-network structure and intervention timing) rather than overstating precision.
  *Notes:* This R0 describes population-level spread within the affected sexual network during the 2022 outbreak, NOT general-population respiratory R0 (which is much lower -- mpox is comparatively inefficient at close-contact/respiratory spread outside of prolonged skin contact). Scenario-level parameters distinguish sexual-contact from non-sexual/household spread explicitly (see sexual_contact_secondary_attack_rate below).
- **Household (non-sexual close-contact) secondary attack rate:** Beta distribution, standard for bounded probability. Point estimate anchored to the systematic-review pooled figure (8%); interval widened to include both the lower UK clade-IIb estimate (4%) and the higher WHO-reported clade-Ib figure (~20%), since clade and setting materially change this number.
  *Notes:* This parameter describes NON-SEXUAL household contact only. See sexual_contact_secondary_attack_rate for the very different sexual-contact figure (~75%), which is the operative parameter for the gathering-associated scenario.
- **Secondary attack rate among sexual contacts (gathering-associated context):** Reported as a single point estimate (~75%) without a formal CI in the source WHO situation report; low/high bounds here are an analyst-assigned +/-15pp uncertainty band pending a primary study with a reported interval, and should be treated as approximate framing rather than a precise statistical CI.
  *Notes:* Directly informs the mpox_gathering scenario's within-event transmission risk.

## Norovirus (GI/GII genogroups)

*enteric virus (fecal-oral, fomite, and short-range aerosol from vomiting). Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Effective reproduction number, person-to-person route, school outbreak setting | 8.92 | [4, 14] | dimensionless | Single study | Xu et al. 2023, Jiangsu Province school outbreaks: effective reproduction number for hu... |
| Incubation period | 1.37 | [0.5, 2] | days | Meta-analysis | Chan et al. meta-analysis of 1022 outbreaks (mean 32.8h [95% CI 30.9-34.6h], median 33.... |
| Symptomatic illness duration | 1.84 | [0.5, 6] | days | Meta-analysis | Chan et al. meta-analysis (mean symptomatic period 44.2h [95% CI 38.9-50.7h]); illness ... |
| Negative-binomial dispersion parameter (k) | 0.5 | [0.15, 1.5] | dimensionless | ⚠️ LOW CONFIDENCE | No dedicated norovirus k-estimation study located in this search pass |
| Basic reproduction number (R0), general community | 2.1 | [1.1, 3.7] | dimensionless | ⚠️ LOW CONFIDENCE | General epidemiological reference range for community norovirus transmission; no single... |
| Outbreak attack rate in closed settings (school/camp proxy for household SAR) | 0.3 | [0.21, 0.42] | probability | Multi-study | Torner et al. 2022 (Catalonia schools/summer camps, 2017-2019): overall attack rate 21.... |

**Justifications and notes:**

- **Effective reproduction number, person-to-person route, school outbreak setting:** Single-study estimate, but directly relevant: this is exactly the "closed indoor setting, sustained person-to-person spread" scenario this project models for norovirus, unlike the more general community R0 above. Low/high bounds are an analyst-assigned uncertainty band (no CI reported in source) pending replication.
  *Notes:* Striking illustration of how much higher effective R can run in closed settings relative to general community R0 -- directly motivates this project's scenario-level contact_rate_multiplier design.
- **Incubation period:** Strongest-evidence parameter in this pathogen's set: large meta-analysis (n=1022 outbreaks) with tight confidence intervals. Values converted from hours to days.
- **Symptomatic illness duration:** Same meta-analysis as incubation period, for internal consistency; converted from hours to days.
  *Notes:* Viral SHEDDING continues far longer than symptomatic illness -- one cited study found shedding starting ~36h post-infection and continuing ~26 days on average (range 11-54 days). This project's infectious_period represents the high-transmission-risk symptomatic window used for the branching-process model; the much longer low-level shedding tail is a known simplification (see docs/validation_plan.md limitations).
- **Negative-binomial dispersion parameter (k):** PLACEHOLDER ASSUMPTION. Norovirus outbreaks are frequently driven by point-source contamination events (a single contaminated food item or surface) that can produce apparent clustering resembling superspreading without true person-to-person transmission heterogeneity -- this confound makes borrowing a k estimate from person-to-person-only pathogens inappropriate, but no norovirus-specific study was found to use instead. Treat this value as a rough placeholder only.
- **Basic reproduction number (R0), general community:** HONESTLY FLAGGED WEAKER PARAMETER: this project could not locate a primary-literature community-baseline R0 estimation study for norovirus in the time available (searches instead surfaced strong data on closed-setting effective R, incubation, and outbreak attack rates -- see closed_setting_effective_r and secondary_attack_rate below, which rest on firmer evidence). Point estimate/range reflects commonly cited textbook figures; treat with caution for any use where the general-community baseline (as opposed to a specific closed-setting outbreak) matters.
  *Notes:* Re-derive from primary literature (e.g. a dedicated norovirus transmission-modeling review) before using this figure for anything consequential.
- **Outbreak attack rate in closed settings (school/camp proxy for household SAR):** NOTE ON DEFINITION: this is a total-outbreak attack rate (cases / total exposed population in a closed setting), not a classic dyadic household secondary attack rate -- no precise per-contact household SAR study was located in this search pass. Used here as the best available proxy since the modeled scenarios are closed-setting outbreaks, not household transmission chains.
  *Notes:* closed_setting_effective_r below provides a complementary, independently-sourced cross-check on transmission intensity in exactly this kind of setting.

## SARS-CoV-2 (ancestral / wild-type strain, 2020 baseline)

*respiratory virus (airborne + droplet transmission). Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Generation interval (time between successive infections) | 5.2 | [3.95, 6.4] | days | Multi-study | Ganyani et al. 2020 Euro Surveill (Singapore cluster: 5.20 days; Tianjin cluster: 3.95 ... |
| Incubation period | 5.8 | [5.01, 6.69] | days | Meta-analysis | Meta-analytic lognormal incubation period (mean 5.8 days, meanlog=1.63, sdlog=0.50), as... |
| Infectious period (duration of transmission risk) | 8 | [5, 10] | days | Expert consensus | Synthesized from generation-interval literature (He et al. 2020 Nat Med; Ganyani et al.... |
| Negative-binomial dispersion parameter (k) | 0.15 | [0.04, 0.6] | dimensionless | Meta-analysis | Endo et al. 2020 Wellcome Open Research (k~0.1); Adam et al. 2020 Nat Med, Hong Kong (k... |
| Basic reproduction number (R0) | 2.66 | [2.41, 2.94] | dimensionless | Meta-analysis | Alimohamadi et al. 2020, meta-analysis of early R0 estimates (pooled R0=2.66, 95% CI 2.... |
| Household secondary attack rate | 0.19 | [0.05, 0.3] | probability | Meta-analysis | Madewell et al. 2020/2022 JAMA Netw Open household-transmission meta-analyses (pooled ~... |

**Justifications and notes:**

- **Generation interval (time between successive infections):** Gamma is the distribution family used directly by the source studies for generation interval; used here to convert branching-process generations into calendar time for epidemic-curve visualization, and as the serial-interval proxy in transmission_inference/epi_reconstruction.py.
  *Notes:* Two independent cluster studies gave meaningfully different point estimates (3.95 vs 5.20 days), reflecting real between-context variation, not just sampling noise.
- **Incubation period:** Lognormal is the distribution family used by essentially all major incubation-period meta-analyses for COVID-19 (right-skewed, long tail consistent with rare long incubations). meanlog/sdlog values are reported directly by source, avoiding re-derivation error.
  *Notes:* Multiple independent meta-analyses converge in the 5.1-6.6 day range for the ancestral strain (Elias et al. 2021: 6.4 days [2.3-17.6]; Wu et al. 2022 JAMA Netw Open: 6.57 days pooled across 141 studies). Omicron incubation is substantially shorter (pooled ~3.4-3.5 days, Wu et al. 2022; Manica et al. 2023 BMC Medicine). Use the ancestral-strain value as default; swap in a variant-specific estimate for contemporary-outbreak scenarios.
- **Infectious period (duration of transmission risk):** No single study reports "infectious period" as a directly measured quantity comparable to incubation period; it is reconstructed from viral-shedding and generation-interval studies. Gamma distribution matches convention used for generation intervals in the cited sources.
  *Notes:* LOWER-CONFIDENCE PARAMETER: distinguish from generation interval (time between successive infections in a chain, commonly cited ~5.2 days for the ancestral strain per Ganyani et al.), which is a related but distinct quantity used directly in transmission-inference calculations (see transmission_inference module). The 8-day figure here represents the outer bound of meaningful transmission risk (pre-symptomatic window + ~5-8 days post-onset for mild cases); actual per-day infectiousness is front-loaded near symptom onset, not uniform across this window -- the branching process model treats this as a generation-interval-driven discrete-generation process rather than a continuous infectious window, which sidesteps needing a precise "period" cutoff (see models/branching_process.py).
- **Negative-binomial dispersion parameter (k):** Lognormal chosen because k estimates in the literature span roughly an order of magnitude and are typically summarized/compared on a log scale. Central estimate anchored to the two most-cited primary estimates (Endo et al., Adam et al.); the [0.04, 0.60] interval reflects the spread found across the 28-study systematic review rather than any single study's narrower CI, since context (setting, time period, level of concurrent control measures) clearly shifts k substantially.
  *Notes:* This is THE key superspreading parameter: k=0.1 implies roughly the "20% of cases cause 80% of transmission" rule of thumb. Setting-specific estimates from the same systematic review ranged 0.014-0.72; New Zealand genomic-cluster-based estimate was notably higher (k=0.63, Tran-Kiem & Bedford), attributed to stronger contact-tracing-era suppression reducing superspreading opportunity. Treat k as context-dependent, not a pathogen constant -- scenario-level overrides are supported (see scenarios.yaml).
- **Basic reproduction number (R0):** Pooled estimate across many early-pandemic studies; used as the community-transmission baseline that scenario-specific contact_rate_multiplier values (see scenarios.yaml) scale up from. The same meta-analysis notes very high heterogeneity between studies (small-study effects, method-dependence) -- treat the interval as indicative, not exact.
  *Notes:* Variant- and context-dependent: later variants (Delta, Omicron) had substantially higher transmissibility; superspreading-driven early growth-phase analyses have estimated R0 as high as 4.7-11.4 in some settings (Sneppen et al. 2021, PMC7540800) once superspreader dynamics are accounted for explicitly. This project treats such elevated values as *scenario-level effective R*, not baseline R0 -- see the contact_rate_multiplier field on each SARS-CoV-2 scenario.
- **Household secondary attack rate:** Beta distribution is the standard choice for a bounded [0,1] probability; parameters implied by the point estimate and CI width via method-of-moments (see models/distributions.py:beta_from_mean_ci).
  *Notes:* Ancestral/Delta-era estimate. Household SAR increased substantially for Omicron (~35% per cross-pathogen benchmarking sources) due to immune escape, and dropped for SARS-CoV-1 (~6%). This project's default scenarios use the ancestral-strain value for consistency with the 2020 case studies (Skagit choir, military basic training) used for external validation.

## Varicella (chickenpox, primary VZV infection)

*respiratory virus (airborne) + direct vesicular contact. Last reviewed: 2026-07-20.*

| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |
|---|---|---|---|---|---|
| Incubation period | 15 | [10, 21] | days | Meta-analysis | CDC Varicella Outbreak Control and Investigation Manual; CDC Pink Book; CDC Yellow Book... |
| Infectious period (relative to rash onset) | 6 | [5, 9] | days | Expert consensus | CDC clinical guidance: infectious 1-2 days before rash onset through full lesion crusti... |
| Negative-binomial dispersion parameter (k) | 1 | [0.3, 2] | dimensionless | ⚠️ LOW CONFIDENCE | No dedicated varicella k-estimation study located in this search pass |
| Basic reproduction number (R0) | 10 | [7, 12] | dimensionless | ⚠️ LOW CONFIDENCE | Commonly cited textbook range for varicella R0 (order of magnitude below measles, above... |
| Household secondary attack rate, susceptible contacts | 0.85 | [0.65, 0.9] | probability | Meta-analysis | CDC Yellow Book 2024 (~85%, range 61-100%); Arvin 1996 / Ross 1962 (90%, as cited in Ma... |

**Justifications and notes:**

- **Incubation period:** Extremely well-replicated figure across independent CDC reference documents spanning many years.
- **Infectious period (relative to rash onset):** Textbook-stable figure, consistent across CDC reference documents.
- **Negative-binomial dispersion parameter (k):** PLACEHOLDER ASSUMPTION by analogy to measles (see measles.k_dispersion.justification for the same reasoning applied there). Not a literature-derived estimate.
  *Notes:* Highest-priority parameter to replace with a real literature estimate if this project is extended.
- **Basic reproduction number (R0):** HONESTLY FLAGGED: unlike measles, no primary modern R0-estimation study for varicella was retrieved in this search pass. Value reflects standard epidemiological teaching (varicella as a highly transmissible but somewhat-less-efficient airborne spreader than measles) rather than a specific citable point estimate. Re-derive from a dedicated source (e.g. Anderson & May-style historical estimates) before using this figure for anything consequential.
- **Household secondary attack rate, susceptible contacts:** Beta distribution for bounded probability; three independent CDC/review sources converge tightly, giving high confidence in this figure specifically for household contacts.
  *Notes:* IMPORTANT COUNTEREXAMPLE TO 'closed settings always increase SAR': a military training centre outbreak investigation (Bhatta et al. 2022) found a much lower SAR of 21.43% among close contacts in that setting, plausibly reflecting partial population immunity (prior infection/vaccination) rather than lower intrinsic transmissibility. This is recorded here as a caution against assuming setting-level multipliers always push SAR upward -- population susceptibility matters at least as much.

## Scenario-level setting parameters

| Scenario | Pathogen | Population | Initial cases | Contact multiplier | Real-world benchmark(s) |
|---|---|---|---|---|---|
| Indoor choir rehearsal (superspreading event) | sars_cov_2 | 61 | 1 | 3.2x | 86.0% |
| Seasonal influenza outbreak in a closed institutional setting | influenza | 300 | 2 | 1.3x | 20.0% |
| Measles outbreak in an under-vaccinated school/community setting | measles | 400 | 1 | 1.0x | 85.0% |
| Military basic training congregate housing | sars_cov_2 | 240 | 2 | 2.0x | 22.3%; 12.8% |
| Mpox transmission at a social/sexual-network-connected gathering | mpox | 40 | 1 | 1.0x | 75.0% |
| Norovirus outbreak in a closed institutional setting (school/camp) | norovirus | 150 | 1 | 1.0x | 21.3%; 33.4% |
| University residence hall outbreak | sars_cov_2 | 600 | 3 | 1.5x | 14.0%; 6.2% |
| Varicella outbreak in a school/childcare setting | varicella | 200 | 1 | 1.0x | 85.0%; 21.4% |
