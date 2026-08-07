"""
Schema definitions for epidemiological parameters and outbreak scenarios.

Design rationale
-----------------
Every quantitative claim in this project (an R0, an incubation period, an
attack rate) is data, not code. Data lives in versioned YAML files
(see data/parameters/ and data/scenarios/) and is loaded through these
Pydantic schemas, which enforce:

  1. Type correctness      (a duration must be a number, not a string)
  2. Range/domain validity  (a probability must be in [0, 1]; a duration
                             must be positive; k must be positive)
  3. Provenance             (every estimate must carry a source and,
                             where applicable, an uncertainty range —
                             an estimate with no citation fails validation)
  4. Referential integrity  (a scenario that references "measles" must
                             have a matching pathogen entry; interventions
                             referenced by name must exist)

This is deliberately stricter than most research code. The cost of a
silently-wrong parameter (say, an R0 typed as 25 instead of 2.5) propagating
through 100,000 Monte Carlo draws is high, and schema validation is cheap
insurance against it.
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field, field_validator, model_validator

# --------------------------------------------------------------------------
# Enumerations
# --------------------------------------------------------------------------

class DistributionFamily(str, Enum):
    """Probability distribution families used throughout the project.

    See docs/design_document.md ("Statistical Framework") for the rationale
    behind each mapping of parameter type -> distribution family.
    """

    BETA = "beta"                    # bounded [0,1]: probabilities, SARs, efficacies
    GAMMA = "gamma"                  # positive, right-skewed: durations, generation intervals
    LOGNORMAL = "lognormal"          # positive, right-skewed: incubation periods
    NEGATIVE_BINOMIAL = "negative_binomial"  # counts with overdispersion: offspring distribution
    POISSON = "poisson"              # counts, homogeneous-mixing limit (k -> infinity)
    UNIFORM = "uniform"              # bounded, non-informative: used only as an explicit fallback
    POINT = "point"                  # degenerate / fixed value: used for deterministic sensitivity runs


class EvidenceQuality(str, Enum):
    """How much weight a parameter estimate should be given.

    This field exists specifically so the project never silently treats a
    single case study's number the same way it treats a pooled
    meta-analytic estimate. See docs/validation_plan.md.
    """

    META_ANALYSIS = "meta_analysis"          # pooled systematic review / meta-analysis
    MULTI_STUDY = "multi_study"              # >=2 independent primary studies, informally pooled
    SINGLE_STUDY = "single_study"            # one primary study or outbreak investigation
    EXPERT_CONSENSUS = "expert_consensus"    # textbook / public-health-agency guidance, no single citable study
    MODELED_ESTIMATE = "modeled_estimate"    # inferred by a transmission model fitted to case counts
    LOW_CONFIDENCE = "low_confidence"        # sparse, conflicting, or indirect evidence; flagged explicitly


# --------------------------------------------------------------------------
# Parameter-level schema
# --------------------------------------------------------------------------

class ParameterEstimate(BaseModel):
    """A single sourced, quantified epidemiological parameter.

    This is the atomic unit of the evidence table required by the project
    spec ("parameter / value / uncertainty range / source / justification").
    """

    name: str = Field(..., description="Machine-readable parameter name, e.g. 'r0'")
    display_name: str = Field(..., description="Human-readable label")
    distribution: DistributionFamily
    point_estimate: float = Field(..., description="Central estimate (mean or median, per source)")
    low: float | None = Field(None, description="Lower bound of reported uncertainty interval")
    high: float | None = Field(None, description="Upper bound of reported uncertainty interval")
    ci_level: float | None = Field(0.95, description="Confidence/credible level of (low, high), e.g. 0.95")
    unit: str = Field(..., description="Unit of measurement, e.g. 'days', 'dimensionless', 'probability'")
    source: str = Field(..., min_length=3, description="Citation: authors, year, journal/agency")
    source_url: str | None = None
    evidence_quality: EvidenceQuality
    justification: str = Field(..., min_length=10, description="Why this value/distribution was chosen")
    justification_detail: str | None = None
    notes: str | None = Field(None, description="Caveats, controversy, context-dependence")

    @field_validator("point_estimate")
    @classmethod
    def _finite(cls, v: float) -> float:
        if v != v or v in (float("inf"), float("-inf")):  # NaN / inf check without numpy dependency
            raise ValueError("point_estimate must be a finite number")
        return v

    @model_validator(mode="after")
    def _bounds_consistent(self) -> ParameterEstimate:
        if self.low is not None and self.high is not None and self.low > self.high:
            raise ValueError(f"{self.name}: low ({self.low}) > high ({self.high})")
        if self.low is not None and self.high is not None:
            if not (self.low <= self.point_estimate <= self.high):
                raise ValueError(
                    f"{self.name}: point_estimate {self.point_estimate} outside "
                    f"reported interval [{self.low}, {self.high}]"
                )
        return self

    @model_validator(mode="after")
    def _domain_by_unit(self) -> ParameterEstimate:
        """Range checks tied to the semantic unit of the parameter."""
        if self.unit == "probability":
            for label, val in [("point_estimate", self.point_estimate), ("low", self.low), ("high", self.high)]:
                if val is not None and not (0.0 <= val <= 1.0):
                    raise ValueError(f"{self.name}.{label}={val} is not a valid probability in [0,1]")
        if self.unit in ("days", "hours"):
            for label, val in [("point_estimate", self.point_estimate), ("low", self.low), ("high", self.high)]:
                if val is not None and val <= 0:
                    raise ValueError(f"{self.name}.{label}={val} must be > 0 for a duration")
        if self.name == "k_dispersion":
            for label, val in [("point_estimate", self.point_estimate), ("low", self.low), ("high", self.high)]:
                if val is not None and val <= 0:
                    raise ValueError(f"{self.name}.{label}={val} (dispersion k) must be > 0")
        return self


class PathogenParameterSet(BaseModel):
    """All sourced parameters for a single pathogen."""

    pathogen_id: str = Field(..., description="Machine-readable id, e.g. 'sars_cov_2'")
    display_name: str
    pathogen_class: str = Field(..., description="e.g. 'respiratory virus', 'enteric virus'")
    parameters: dict[str, ParameterEstimate]
    last_reviewed: str = Field(..., description="ISO date this parameter set was last checked against literature")
    review_method: str = Field(
        ...,
        description=(
            "How this evidence was gathered. Must state explicitly if this is a targeted "
            "literature search rather than a PRISMA-compliant systematic review — "
            "see docs/evidence_tables.md for the methodology statement."
        ),
    )

    @model_validator(mode="after")
    def _required_core_parameters(self) -> PathogenParameterSet:
        required = {"r0", "incubation_period", "infectious_period", "k_dispersion", "secondary_attack_rate"}
        missing = required - set(self.parameters.keys())
        if missing:
            raise ValueError(f"{self.pathogen_id}: missing required core parameters {missing}")
        return self


# --------------------------------------------------------------------------
# Scenario-level schema
# --------------------------------------------------------------------------

class PopulationStructure(BaseModel):
    """Describes the contact/population structure of an outbreak setting."""

    population_size: int = Field(..., gt=0)
    initial_cases: int = Field(1, ge=1)
    setting_type: str = Field(..., description="e.g. 'indoor gathering', 'congregate housing'")
    contact_rate_multiplier: float = Field(
        1.0, gt=0,
        description=(
            "Multiplier applied to the pathogen's baseline R0 to reflect this setting's "
            "elevated (or reduced) contact intensity, duration, and environmental risk "
            "relative to average community transmission. 1.0 = community baseline."
        ),
    )
    contact_rate_multiplier_low: float | None = Field(None, gt=0)
    contact_rate_multiplier_high: float | None = Field(None, gt=0)
    exposure_duration_hours: float | None = Field(None, gt=0)
    ventilation_ach: float | None = Field(None, gt=0, description="Baseline air changes per hour, if known/assumed")


class ObservedOutcome(BaseModel):
    """A real, literature-reported outcome used for external validation."""

    description: str
    attack_rate: float | None = Field(None, ge=0, le=1)
    outbreak_size: int | None = Field(None, ge=0)
    population_at_risk: int | None = Field(None, gt=0)
    source: str
    source_url: str | None = None
    caveats: str | None = Field(
        None, description="Known controversies or limitations of this specific data point"
    )


class OutbreakScenario(BaseModel):
    """A named, reproducible outbreak configuration: pathogen + setting + (optional) real-world benchmark."""

    scenario_id: str
    display_name: str
    pathogen_id: str
    description: str
    population: PopulationStructure
    observed_outcomes: list[ObservedOutcome] = Field(default_factory=list)
    assumptions: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)

    @field_validator("scenario_id")
    @classmethod
    def _id_format(cls, v: str) -> str:
        if not v.replace("_", "").isalnum():
            raise ValueError("scenario_id must be alphanumeric/underscore")
        return v
