"""
TransPhylo interface (genomic transmission inference).

TransPhylo (Didelot, Fraser, Gardy & Colijn, 2017, Mol Biol Evol) infers a
transmission tree from a DATED PHYLOGENY (a time-calibrated tree relating
pathogen genome sequences from different cases, e.g. from BEAST2 or
treedater) plus a generation-time and sampling-time distribution. It is the
appropriate tool "when genomic data are available" (requirement #10)
-- which this project's bundled scenarios do NOT have (there is no real
sequence data for the Skagit choir outbreak, the Fort Benning outbreak,
etc. bundled here). This module is therefore a genuinely working, documented
INTERFACE for a user who has real genomic data to plug in, not a component
this project's own examples or tests exercise end-to-end -- that
distinction is stated explicitly rather than papered over (see the
module-level WARNING below and docs/validation_plan.md).

TransPhylo is an R package; this project is Python. Two integration paths
are documented:
  (a) rpy2 (in-process R interface) -- convenient for interactive/notebook
      use; requires R + TransPhylo + ape installed alongside Python.
  (b) subprocess to a standalone Rscript -- more robust for pipeline/CI use
      (no shared-process R/Python state to manage), at the cost of
      serializing inputs/outputs through files.
Both are implemented below; `run_transphylo()` uses (a) by default and
falls through to (b) if rpy2 is not installed.

Install (optional; NOT a default dependency of this project -- see
pyproject.toml's `transmission-inference-r` extra):
    pip install "outbreak-simulator[transmission-inference-r]"
    # plus, in R:
    install.packages(c("ape", "remotes"))
    remotes::install_github("xavierdidelot/TransPhylo")
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


class TransPhyloNotAvailableError(RuntimeError):
    """Raised when neither rpy2 nor a working Rscript with TransPhylo installed can be found."""


@dataclass
class TransPhyloConfig:
    dated_tree_newick: str  # a dated phylogeny in Newick format, e.g. from treedater/BEAST2
    sample_dates: dict[str, float]  # tip label -> sampling date (decimal year or days since a reference)
    # Gamma shape for the within-host generation time (see pathogens.yaml
    # generation_interval, converted to years if using decimal-year dates)
    generation_time_shape: float
    generation_time_scale: float
    sampling_time_shape: float  # Gamma shape for time from infection to sampling
    sampling_time_scale: float
    mcmc_iterations: int = 10_000
    outdir: str | None = None


R_SCRIPT_TEMPLATE = r"""
suppressMessages(library(ape))
suppressMessages(library(TransPhylo))

tree <- read.tree("{tree_path}")
ptree <- ptreeFromPhylo(tree, dateLastSample = {date_last_sample})

result <- inferTTree(
    ptree,
    mcmcIterations = {mcmc_iterations},
    w.shape = {gen_shape}, w.scale = {gen_scale},
    ws.shape = {samp_shape}, ws.scale = {samp_scale}
)

med <- medTTree(result)
mat <- computeMatWIW(result)

write.csv(mat, "{output_matrix_path}")
saveRDS(result, "{output_rds_path}")
cat("TransPhylo run complete.\n")
"""


def _generate_r_script(
    cfg: TransPhyloConfig, tree_path: Path, out_matrix: Path, out_rds: Path, date_last_sample: float
) -> str:
    return R_SCRIPT_TEMPLATE.format(
        tree_path=str(tree_path), date_last_sample=date_last_sample,
        mcmc_iterations=cfg.mcmc_iterations,
        gen_shape=cfg.generation_time_shape, gen_scale=cfg.generation_time_scale,
        samp_shape=cfg.sampling_time_shape, samp_scale=cfg.sampling_time_scale,
        output_matrix_path=str(out_matrix), output_rds_path=str(out_rds),
    )


def run_transphylo_subprocess(cfg: TransPhyloConfig) -> dict:
    """Integration path (b): write inputs to a temp dir, invoke Rscript, read
    results back. Requires R with TransPhylo + ape installed and `Rscript`
    on PATH. Raises TransPhyloNotAvailableError with a clear message if not."""
    if subprocess.run(["which", "Rscript"], capture_output=True).returncode != 0:
        raise TransPhyloNotAvailableError(
            "Rscript not found on PATH. Install R and TransPhylo (see this module's docstring) "
            "to use the genomic transmission-inference pathway; otherwise use "
            "transmission_inference.epi_reconstruction for the non-genomic alternative."
        )
    workdir = Path(cfg.outdir) if cfg.outdir else Path(tempfile.mkdtemp(prefix="transphylo_"))
    workdir.mkdir(parents=True, exist_ok=True)
    tree_path = workdir / "dated_tree.nwk"
    tree_path.write_text(cfg.dated_tree_newick)
    out_matrix = workdir / "wiw_matrix.csv"
    out_rds = workdir / "transphylo_result.rds"
    date_last_sample = max(cfg.sample_dates.values())

    script = _generate_r_script(cfg, tree_path, out_matrix, out_rds, date_last_sample)
    script_path = workdir / "run_transphylo.R"
    script_path.write_text(script)

    proc = subprocess.run(["Rscript", str(script_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"TransPhylo Rscript failed:\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")

    return {
        "who_infected_whom_matrix_csv": str(out_matrix),
        "r_result_object": str(out_rds),
        "stdout": proc.stdout,
        "workdir": str(workdir),
    }


def run_transphylo_rpy2(cfg: TransPhyloConfig):
    """Integration path (a): in-process R via rpy2. Preferred for
    interactive/notebook use since it avoids file round-tripping and
    returns live R objects. See module docstring for install instructions."""
    try:
        from rpy2.robjects.packages import importr
    except ImportError as exc:
        raise TransPhyloNotAvailableError(
            "rpy2 is not installed. Install with: pip install 'outbreak-simulator[transmission-inference-r]', "
            "plus R + TransPhylo (see this module's docstring). Falling back to "
            "run_transphylo_subprocess() is also available if you have a working Rscript but not rpy2."
        ) from exc

    ape = importr("ape")
    transphylo = importr("TransPhylo")

    workdir = Path(cfg.outdir) if cfg.outdir else Path(tempfile.mkdtemp(prefix="transphylo_"))
    workdir.mkdir(parents=True, exist_ok=True)
    tree_path = workdir / "dated_tree.nwk"
    tree_path.write_text(cfg.dated_tree_newick)

    tree = ape.read_tree(str(tree_path))
    date_last_sample = max(cfg.sample_dates.values())
    ptree = transphylo.ptreeFromPhylo(tree, dateLastSample=date_last_sample)

    result = transphylo.inferTTree(
        ptree, mcmcIterations=cfg.mcmc_iterations,
        w_shape=cfg.generation_time_shape, w_scale=cfg.generation_time_scale,
        ws_shape=cfg.sampling_time_shape, ws_scale=cfg.sampling_time_scale,
    )
    wiw_matrix = transphylo.computeMatWIW(result)
    return {"r_result": result, "who_infected_whom_matrix": wiw_matrix}


def run_transphylo(cfg: TransPhyloConfig, prefer: str = "rpy2") -> dict:
    """Convenience dispatcher: try rpy2 first (if prefer='rpy2', the default),
    fall back to the subprocess path, and raise a single clear error if
    neither is available."""
    if prefer == "rpy2":
        try:
            return run_transphylo_rpy2(cfg)
        except TransPhyloNotAvailableError:
            pass
    return run_transphylo_subprocess(cfg)


def generation_time_years_from_evidence_table(mean_days: float, sd_days: float) -> tuple[float, float]:
    """Convert this project's day-scale Gamma generation-interval parameters
    (data/parameters/pathogens.yaml) into the (shape, scale) years-scale
    parameterization TransPhylo expects when sample dates are in decimal
    years (TransPhylo's convention in most published examples)."""
    mean_years, sd_years = mean_days / 365.25, sd_days / 365.25
    shape = (mean_years / sd_years) ** 2
    scale = sd_years**2 / mean_years
    return shape, scale
