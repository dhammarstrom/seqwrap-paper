##############################################################################
#
# Make documentation
#
# This file runs the full analysis pipeline (model fitting, figures) and
# renders the manuscript and supplement. If model estimates are not present in
# analysis/data/derived_data/, functions run_modelN will place reproduced
# estimates there. Rerunning all model estimates will take > 2 hours.
# To overwrite previous results, set overwrite_models to TRUE (row 70).
#
# Simulation results are re-run if make_sim is set to TRUE. Running simulation
# results takes > 24 hours (dependent on number of cores). We recommend
# downloading simulation results from Dataverse
# (https://doi.org/10.18710/I7U71O), see the README for instructions.
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
# make_sim          re-run both simulation scenarios. > 24 hours.
# make_meth_perm    re-run the methylation permutation study. Hours.
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
make_meth_perm <- FALSE

# Detect cores
cores <- parallel::detectCores()


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
} else {
  download_dataverse()
}

# Check if simulation results are present
if (!dir.exists(here::here("analysis/data/raw_data"))) {
  warning(
    "Simulation results are not present in this repository. ",
    "Results can be downloaded from Dataverse ",
    "(https://doi.org/10.18710/I7U71O). The simulation results are needed ",
    "to reproduce results in the manuscript."
  )
}



# Methylation permutation study ----------------------------------------------
#
# analysis/R/methylation-error-permutation.R estimates type I error rates for
# the beta and M-value models from within-participant permuted ("plasmode")
# data: 500 permutations x 1000 sites x 4 model arms, written one file per
# permutation into analysis/data/derived_data/permutation_v2/.
#
# The script skips a permutation whose file already exists, so an interrupted
# job resumes simply by running it again. Note that permutation_v2/ is a NEW
# directory: permutation/ holds the earlier two-arm run, kept because its
# evaluation columns use the old names and the two cannot be pooled.
#
# case-study-methylation.qmd sources this script itself when the directory is
# empty. Running it here instead makes a job of several hours a deliberate
# step rather than a side effect of rendering a document.
#
# Inputs are produced by analysis/R/methylation-case-study-dataprep.R, which
# downloads the raw arrays from GEO and writes the rgset, the sample metadata
# and the normalized, probe-filtered gset.

if (make_meth_perm) {
  der <- here::here("analysis/data/derived_data")
  meth_inputs <- file.path(
    der,
    c("seaborne-gset-normalized.RDS", "seaborne-metadata.RDS")
  )

  # Data preparation. Every DOWNLOAD in that script is guarded by
  # file.exists(): the SOFT series metadata (pheno_raw.rds), the ~757 MB
  # GSE114763_RAW.tar, and the untar/gunzip of the IDATs. So is the rgset
  # build and the normalization. Re-running never re-downloads and never
  # redoes work that is already on disk; the source() below is skipped
  # entirely once both outputs exist.
  if (!all(file.exists(meth_inputs))) {
    source(here::here("analysis/R/methylation-case-study-dataprep.R"))
  }

  if (!all(file.exists(meth_inputs))) {
    stop(
      "Data preparation did not produce: ",
      paste(basename(meth_inputs[!file.exists(meth_inputs)]), collapse = ", ")
    )
  }

  source(here::here("analysis/R/methylation-error-permutation.R"))
}


# Source figure files (these are needed for the manuscript and supplement)
source(here::here("analysis/figures/figure-2.R"))
source(here::here("analysis/figures/figure-3.R"))
source(here::here("analysis/figures/figure-4.R"))

# Render documentation
quarto::quarto_render(here::here("analysis/paper/supplement.qmd"))
quarto::quarto_render(here::here("analysis/paper/paper.qmd"))
