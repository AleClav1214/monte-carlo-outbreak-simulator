"""
Distribution helpers.

Bridges the epidemiological parameterization used in the evidence tables
(point estimate + uncertainty interval, in domain-meaningful units) to the
scipy/numpy objects the simulation engine actually samples from.

Why each distribution family is used where (see also docs/design_document.md,
"Statistical Framework"):

  Beta        - any bounded [0,1] probability: secondary attack rates,
                intervention effectiveness, compliance rates, vaccine
                efficacy. Natural bounded support; two-parameter family
                fits a mean + a credible interval width without extra
                degrees of freedom.
  Gamma       - positive, right-skewed durations: infectious periods,
                generation/serial intervals. This is the convention used
                directly by the primary literature we cite (e.g. Ganyani
                et al. report generation intervals as Gamma-distributed).
  Lognormal   - positive, right-skewed durations with a heavier tail than
                Gamma is typically fit with: incubation periods. Again
                matches literature convention directly (e.g. Lauer et al.
                report incubation period as Lognormal with meanlog/sdlog).
  NegBinomial - individual-level offspring (secondary case) counts. This
                is THE standard mechanism for representing superspreading
                (Lloyd-Smith et al. 2005): mean = R, dispersion = k, and
                k -> infinity recovers the Poisson (homogeneous-mixing)
                limit.
  Poisson     - homogeneous-mixing count processes, used inside the
                compartmental (SEIR) model's stochastic transition counts,
                where individual-level heterogeneity is not represented.

Every function below samples parameter UNCERTAINTY (the outer Monte Carlo
loop: "what is the true population R0, given what the literature says?"),
which is distinct from the individual-level STOCHASTICITY sampled inside a
single simulation run (the inner loop: "given this run's R0 and k, how many
people does this specific infected person infect?"). Keeping these two
sampling layers conceptually separate is one of this project's central
design decisions -- see docs/design_document.md.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from outbreak_simulator.data.schemas import DistributionFamily, ParameterEstimate


def beta_from_mean_ci(mean: float, low: float, high: float, ci_level: float = 0.95) -> stats.rv_continuous:
    """Fit a Beta distribution's (alpha, beta) via method-of-moments from a
    mean and an approximate CI, using the CI width to imply a variance.

    Method: treat (high - low) as approximately 2 * z * sd for a
    normal-ish approximation at the given ci_level (z=1.96 for 95%), solve
    for sd, then use standard Beta method-of-moments:
        alpha = mean * (mean*(1-mean)/var - 1)
        beta  = (1-mean) * (mean*(1-mean)/var - 1)
    This is an approximation (Beta is not normal), adequate for turning a
    literature-reported interval into a sensible sampling distribution
    without requiring the original study's raw data.
    """
    if not (0 < mean < 1):
        raise ValueError(f"mean must be in (0,1) for a Beta distribution, got {mean}")
    z = stats.norm.ppf(0.5 + ci_level / 2)
    sd = max((high - low) / (2 * z), 1e-6)
    var = sd**2
    # cap variance below the Beta's maximum possible variance mean*(1-mean)
    max_var = mean * (1 - mean) * 0.98
    var = min(var, max_var)
    nu = mean * (1 - mean) / var - 1
    alpha = max(mean * nu, 1e-3)
    beta_param = max((1 - mean) * nu, 1e-3)
    return stats.beta(alpha, beta_param)


def gamma_from_mean_ci(mean: float, low: float, high: float, ci_level: float = 0.95) -> stats.rv_continuous:
    """Fit a Gamma distribution's (shape, scale) via method-of-moments."""
    z = stats.norm.ppf(0.5 + ci_level / 2)
    sd = max((high - low) / (2 * z), 1e-6)
    var = sd**2
    shape = mean**2 / var
    scale = var / mean
    return stats.gamma(a=shape, scale=scale)


def lognormal_from_mean_ci(mean: float, low: float, high: float, ci_level: float = 0.95) -> stats.rv_continuous:
    """Fit a Lognormal distribution's (mu, sigma) via method-of-moments on the
    natural (not log) scale mean/variance, then convert to the log-scale
    parameters scipy expects."""
    z = stats.norm.ppf(0.5 + ci_level / 2)
    sd = max((high - low) / (2 * z), 1e-6)
    var = sd**2
    sigma2 = np.log(1 + var / mean**2)
    sigma = np.sqrt(sigma2)
    mu = np.log(mean) - sigma2 / 2
    return stats.lognorm(s=sigma, scale=np.exp(mu))


def distribution_for_parameter(param: ParameterEstimate) -> stats.rv_continuous:
    """Construct the scipy distribution object a ParameterEstimate implies,
    for sampling parameter-uncertainty draws in the outer Monte Carlo loop."""
    low = param.low if param.low is not None else param.point_estimate * 0.7
    high = param.high if param.high is not None else param.point_estimate * 1.3
    ci = param.ci_level or 0.95

    if param.distribution == DistributionFamily.BETA:
        return beta_from_mean_ci(param.point_estimate, low, high, ci)
    if param.distribution == DistributionFamily.GAMMA:
        return gamma_from_mean_ci(param.point_estimate, low, high, ci)
    if param.distribution == DistributionFamily.LOGNORMAL:
        return lognormal_from_mean_ci(param.point_estimate, low, high, ci)
    if param.distribution == DistributionFamily.UNIFORM:
        return stats.uniform(loc=low, scale=(high - low))
    if param.distribution == DistributionFamily.POINT:
        return stats.uniform(loc=param.point_estimate, scale=0.0)
    raise ValueError(
        f"'{param.distribution}' has no direct outer-loop sampling distribution "
        f"(negative_binomial/poisson are individual-level offspring distributions, "
        f"sampled inside a run via negative_binomial_offspring() below, not across runs)"
    )


def sample_parameter(param: ParameterEstimate, rng: np.random.Generator, size: int | None = None):
    """Draw one (or `size`) Monte Carlo sample(s) of a parameter's true value."""
    dist = distribution_for_parameter(param)
    return dist.rvs(size=size, random_state=rng)


def negative_binomial_offspring(
    r_effective: float | np.ndarray, k: float, rng: np.random.Generator, size: int | None = None
) -> np.ndarray:
    """Sample individual-level secondary-case counts from a Negative Binomial
    offspring distribution parameterized by (mean=r_effective, dispersion=k),
    following the standard Lloyd-Smith et al. (2005) parameterization.

    numpy's negative_binomial(n, p) has mean = n(1-p)/p. Solving for the
    (mean=R, dispersion=k) parameterization used throughout the epi
    literature: n=k, p=k/(k+R).
    """
    if k <= 0:
        raise ValueError("k (dispersion) must be > 0")
    r_arr = np.asarray(r_effective, dtype=float)
    if np.any(r_arr < 0):
        raise ValueError("r_effective must be >= 0")
    p = k / (k + r_arr)
    p = np.clip(p, 1e-12, 1.0)  # guard against r_effective == 0 edge case
    return rng.negative_binomial(k, p, size=size)


def poisson_count(rate: float | np.ndarray, rng: np.random.Generator, size: int | None = None) -> np.ndarray:
    """Sample from a Poisson process -- used for the SEIR compartmental
    model's transition counts (homogeneous-mixing limit, k -> infinity)."""
    rate_arr = np.clip(np.asarray(rate, dtype=float), 0, None)
    return rng.poisson(rate_arr, size=size)
