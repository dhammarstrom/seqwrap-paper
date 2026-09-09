##############################################################################
#
# Make documentation
#
# This file runs the full analysis pipeline (model fitting, figures) and
# renders the manuscript and supplement. If model estimates are not present in
# analysis/data/derived_data/, functions run_modelN will place reproduced
# estimates there. Rerunning all model estimates will take > 2 hours.
# To overwrite previous results, set overwrite_models to TRUE (see the
# switches below).
#
# Simulation results are re-run if make_sim is set to TRUE. Running simulation
# results takes > 24 hours (dependent on number of cores). The methylation
# case study is re-run if make_meth is set to TRUE, which takes about two days
# on 16 cores. We recommend downloading both from Dataverse
# (https://doi.org/10.18710/I7U71O), see the README for instructions. When
# either switch is FALSE, download_dataverse() fetches the whole data set and
# places simulation results in analysis/data/raw_data/ and methylation
# results in analysis/data/derived_data/. Files already on disk are skipped.
#
##############################################################################

# Restore the package environment
renv::restore()

# Restore package versions (lmerSeq is locked as a Local source pointing at
# inst/lmerSeq_0.1.7.tar.gz. If installing lmerSeq is problematic it can be
# installed directly from source:
# renv::install("inst/lmerSeq_0.1.7.tar.gz")

# Install the research compendium package
renv::install(".")

# Packages assumed in analysis/R/
library(seqwrappaper)
library(seqwrap)
library(tidyverse)
library(glmmTMB)
library(edgeR)
library(DESeq2)
library(lmerTest)
library(DHARMa)

# Source paper-specific orchestrators (not part of the package API)
source(here::here("analysis/R/model-functions.R"))
source(here::here("analysis/R/simulation-functions.R"))

# ---------------------------------------------------------------------------
# Switches for a full re-run
#
# overwrite_models  refit the five real-data (Pillon) models. > 2 hours.
# make_sim          re-run both simulation scenarios. > 24 hours. When FALSE
#                   the results are downloaded from Dataverse instead.
# make_meth         re-run the methylation case study (data prep, permutation
#                   study, full run). About two days on 16 cores. When FALSE
#                   the results are downloaded from Dataverse instead. Note
#                   that the timing experiment behind Figure 5 F-G
#                   (timing_v1/) and the superseded four-arm run timing
#                   (full_v2/timing.RDS) that the paper quotes are NOT
#                   reproduced by this script under either setting; they are
#                   only available from Dataverse (see below).
#
# overwrite_models MUST be TRUE when re-running against seqwrap >= 0.8.0 with
# cached results produced by an earlier version. The seqwrapResults class
# gained the properties `models`, `targets`, `cache` and `elapsed_time`, and
# S7 validates a deserialised object against the class definition stored
# inside it, so seqwrap_summarise() on an older m1_results.RDS fails with
# "Can't find property <seqwrap::seqwrap_results>@cache". Refitting rewrites
# them in the current format.
# ---------------------------------------------------------------------------

overwrite_models <- FALSE
make_sim <- FALSE
make_meth <- FALSE

# Detect cores
cores <- parallel::detectCores()


# Download precomputed results from Dataverse ---------------------------------
#
# The Dataverse data set (https://doi.org/10.18710/I7U71O) holds both the
# simulation results (placed in analysis/data/raw_data/) and the methylation
# case-study results (labelled derived_data/... on Dataverse and placed in
# analysis/data/derived_data/). download_dataverse() skips files that already
# exist, so this is cheap on a populated checkout.
#
# It runs unconditionally: even with make_sim and make_meth both TRUE, the
# timing experiment (derived_data/timing_v1/) and the superseded four-arm run
# timing (derived_data/full_v2/timing.RDS) are only available from Dataverse,
# and Figure 5 F-G and the paper read them. The re-run steps below overwrite
# whatever else was downloaded.

download_dataverse()


# Fit models on the real-world (Pillon) data ---------------------------------
#
# These run FIRST because the simulation depends on them: the mean-dispersion
# trend models that parameterise simcounts2() are loess fits to model 1's
# dispersion estimates, built in analysis/figures/figure-2.R. Sourcing that
# script before run_model1() has put `m1_results` in the global environment
# fails -- its own fallback is guarded on exists("run_model1"), which is
# always TRUE here because model-functions.R is sourced above, so the fallback
# never fires and figure-2.R aborts on a missing `m1_results`.

m1_results <- run_model1(CORES = cores, overwrite = overwrite_models)
m2_results <- run_model2(CORES = cores, overwrite = overwrite_models)
m3_results <- run_model3(CORES = cores, overwrite = overwrite_models)
m4_results <- run_model4(CORES = cores, overwrite = overwrite_models)
m5_results <- run_model5(CORES = cores, overwrite = overwrite_models)


# Simulations ----------------------------------------------------------------
#
# figure-2.R constructs `trend_model_observed` (inverse-variance weighted, the
# high-variability scenario) and `trend_model_observed_low` (unweighted, the
# low-variability scenario) from m1_results. sim_wrap1 uses the first,
# sim_wrap2 the second.
#
# `overwrite` is passed through deliberately. Both wrappers skip a dataset
# whose results are already on disk, so a re-run over a populated
# analysis/data/raw_data would otherwise do nothing.

if (make_sim) {
  source(here::here("analysis/figures/figure-2.R"))
  sim_wrap1(cores = cores, overwrite = TRUE)
  sim_wrap2(cores = cores, overwrite = TRUE)
}

# Check if simulation results are present
sim_dirs <- here::here(
  "analysis/data/raw_data",
  c("estimates", "estimates2", "evaluations", "evaluations2",
    "simdata", "simdata2")
)
if (!all(dir.exists(sim_dirs))) {
  warning(
    "Simulation results are not present in this repository (missing: ",
    paste(basename(sim_dirs[!dir.exists(sim_dirs)]), collapse = ", "),
    "). Results can be downloaded from Dataverse ",
    "(https://doi.org/10.18710/I7U71O). The simulation results are needed ",
    "to reproduce results in the manuscript."
  )
}


# Methylation case study -----------------------------------------------------
#
# Three scripts, run in this order. Each skips work whose output already
# exists, so an interrupted run resumes by setting make_meth <- TRUE again.
#
# 1. analysis/R/methylation-case-study-dataprep.R
#    Downloads the raw arrays (GSE114763, ~757 MB) from GEO into
#    analysis/data/raw_data/GSE114763/ and writes seaborne-rgset.RDS,
#    seaborne-metadata.RDS and the quantile normalized, probe-filtered
#    seaborne-gset-quantile.RDS to analysis/data/derived_data/. ~15 min.
#    paper.qmd reads pheno_raw.rds and the gset directly.
#
# 2. analysis/R/methylation-error-permutation.R
#    Type I error and power from within-participant permuted data: 200
#    iterations x 4 model arms, one file per iteration, written to
#    analysis/data/derived_data/permutation_v4/. ~25 h on 16 cores.
#    paper.qmd and figure-5.R source it themselves when the directory
#    is empty; running it here makes the job a deliberate step instead.
#
# 3. analysis/R/methylation-case-full.R
#    The two REML arms on all 770 441 sites, written to
#    analysis/data/derived_data/full_v3/ (beta_reml.RDS, m_reml.RDS, their
#    -meta files and timing.RDS). ~20 h on 16 cores. paper.qmd and
#    figure-5.R STOP if the arm files are missing rather than refit.
#
# The two long scripts run in their own R process: the full run holds a
# results object of the order of a gigabyte, and a fresh process per step
# returns memory between steps.
#
# NOT run here: analysis/R/methylation-timing-experiment.R, which produced
# analysis/data/derived_data/timing_v1/ (Figure 5 F-G). It was run against the
# earlier, functionally normalized gset (seaborne-gset-normalized.RDS) and
# validates its extrapolation against full_v2/timing.RDS, the four-arm run on
# those same arrays. Both files are kept in derived_data for that reason; the
# full_v2 arm files are in analysis/archive/. paper.qmd also reads
# full_v2/timing.RDS for the quoted full-run hours so that the text and Figure
# 5 G describe the same run. Timing does not depend on the normalization, so
# the experiment was not repeated on the quantile gset.

if (make_meth) {
  der <- here::here("analysis/data/derived_data")
  meth_inputs <- file.path(
    der,
    c("seaborne-gset-quantile.RDS", "seaborne-metadata.RDS")
  )

  if (!all(file.exists(meth_inputs))) {
    source(here::here("analysis/R/methylation-case-study-dataprep.R"))
  }

  if (!all(file.exists(meth_inputs))) {
    stop(
      "Data preparation did not produce: ",
      paste(basename(meth_inputs[!file.exists(meth_inputs)]), collapse = ", ")
    )
  }

  rscript <- file.path(R.home("bin"), "Rscript")
  for (s in c(
    "analysis/R/methylation-error-permutation.R",
    "analysis/R/methylation-case-full.R"
  )) {
    status <- system2(rscript, here::here(s))
    if (status != 0) stop(s, " failed with status ", status, call. = FALSE)
  }
}

# Check that the methylation case-study inputs are present. paper.qmd and
# figure-5.R read all of these; the array set and metadata unconditionally,
# the full-run arms and timing files with a stop() rather than a refit. The
# permutation study is the one piece that is refitted on the fly when its
# directory is empty (~25 h), so it is checked here as well to make a missing
# download visible before the render starts.
der <- here::here("analysis/data/derived_data")
meth_files <- file.path(
  der,
  c(
    "seaborne-gset-quantile.RDS",
    "seaborne-metadata.RDS",
    "full_v3/beta_reml.RDS",
    "full_v3/m_reml.RDS",
    "full_v3/timing.RDS",
    "full_v2/timing.RDS",
    "timing_v1/timing-cells.RDS",
    "timing_v1/timing-scaling.RDS"
  )
)
meth_missing <- meth_files[!file.exists(meth_files)]
n_perm_files <- length(list.files(
  file.path(der, "permutation_v4"),
  pattern = "^perm_\\d+\\.RDS$"
))
if (length(meth_missing) > 0 || n_perm_files == 0) {
  warning(
    "Methylation case-study results are not present in this repository ",
    "(missing: ",
    paste(
      c(basename(meth_missing), if (n_perm_files == 0) "permutation_v4/"),
      collapse = ", "
    ),
    "). They are downloaded from Dataverse (https://doi.org/10.18710/I7U71O) ",
    "by download_dataverse() above, or reproduced with make_meth <- TRUE ",
    "(except timing_v1/ and full_v2/timing.RDS, which are Dataverse only). ",
    "paper.qmd and figure-5.R will fail without them."
  )
}

# Source figure files (these are needed for the manuscript and supplement)
source(here::here("analysis/figures/figure-2.R"))
source(here::here("analysis/figures/figure-3.R"))
source(here::here("analysis/figures/figure-4.R"))
source(here::here("analysis/figures/figure-5.R"))

# Render documentation
quarto::quarto_render(here::here("analysis/paper/supplement.qmd"))
quarto::quarto_render(here::here("analysis/paper/paper.qmd"))
