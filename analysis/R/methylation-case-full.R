# This script fits the full data set from Seaborne et al. using both the beta
# and m models. Raw data is downloaded in methylation-case-study-datapre.R
# Here we proceed from the rgset

## Packages and settings ###################################################
library(GEOquery) # For accessing the raw data
library(minfi) # For pre-processing of raw data
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(maxprobes) # For filtering cross-reactive probes
library(seqwrap)


# Load color scale
source(here::here("analysis/figures/figure-opts.R"))


# Detect number of cores
cores <- parallel::detectCores()


# Quality control and normalization of data ##################################
#
# Moved to analysis/R/methylation-case-study-dataprep.R. The normalized gset is
# an input to both this script and the permutation study
# (analysis/R/methylation-error-permutation.R), so it is built once, upstream,
# rather than as a side effect of the script that fits the full models.

# Load normalized saved gset ################################################

gset_file <- here::here(
  "analysis/data/derived_data/seaborne-gset-normalized.RDS"
)

if (!file.exists(gset_file)) {
  stop(
    "seaborne-gset-normalized.RDS is missing. Run ",
    "analysis/R/methylation-case-study-dataprep.R first -- it downloads the ",
    "arrays if needed and writes the normalized, filtered gset.",
    call. = FALSE
  )
}

gset <- readRDS(gset_file)
metadata <- readRDS(
  here::here("analysis/data/derived_data/seaborne-metadata.RDS")
)

# Deriving beta- and m-values ##############################################

# Calculate the beta-values with the Illumina-style offset (c = 100)
beta_vals <- getBeta(gset, offset = 100)
# Derive M-values from the calculated beta values.
m_vals <- log2(beta_vals / (1 - beta_vals))


# No need to squeeze beta values away from the boundaries: (y(n - 1) + 0.5) / n
# n_obs     <- ncol(beta_vals)
# beta_vals <- (beta_vals * (n_obs - 1) + 0.5) / n_obs

# Clean up
# Remove gset as no longer used
rm(gset)
gc()


# Save data insuitable data frames #########################################
betadf <- data.frame(beta_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

mdf <- data.frame(m_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

# Evaluation function for the full setup ##################################
# This function will give convergence flags. m$fit$convergence flags if
# the numerical procedure is OK (calculation of the log-likelihood).
# see e.g., https://stackoverflow.com/questions/79110546/glmmtmb-convergence-messages
# The pdHess gives TRUE when the the covariance is positive-definite.
# see the above, and
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html

eval_fun <- function(m) {
  data.frame(
    convergence = m$fit$convergence,
    pdHess = m$sdr$pdHess,
    n = nrow(m$frame)
  )
}


# Fit models ################################################################

bm1 <- seqwrap_compose(
  modelfun = glmmTMB::glmmTMB,
  data = betadf,
  metadata = metadata,
  samplename = "geo_accession",
  eval_fun = eval_fun,
  arguments = list(
    formula = y ~ time + (1 | participant),
    family = glmmTMB::beta_family()
  )
)


# Beta-values  model ########################################################

# Both models are cached...
# Check if models results exists. Remove this file if re-run needed.

if (
  !file.exists(
    here::here("analysis/data/derived_data/beta-model-full.RDS")
  )
) {
  bm1_results <- seqwrap(
    bm1,
    return_models = FALSE,
    cores = cores,
    #   subset = 1:10,
    verbose = FALSE,
  )

  saveRDS(
    bm1_results,
    here::here("analysis/data/derived_data/beta-model-full.RDS")
  )

  # Clean-up ##
  rm(bm1_results)
  rm(bm1)
  gc()
}


mm1 <- seqwrap_compose(
  modelfun = glmmTMB::glmmTMB,
  data = mdf,
  metadata = metadata,
  samplename = "geo_accession",
  eval_fun = eval_fun,
  arguments = list(
    formula = y ~ time + (1 | participant)
  )
)


# M-values model ########################################################
if (
  !file.exists(
    here::here("analysis/data/derived_data/m-model-full.RDS")
  )
) {
  mm1_results <- seqwrap(
    mm1,
    return_models = FALSE,
    cores = cores,
    # subset = 1:12,
    verbose = FALSE,
  )

  saveRDS(
    mm1_results,
    here::here("analysis/data/derived_data/m-model-full.RDS")
  )

  # Clean-up ##
  rm(mm1_results)
  rm(mm1)
  gc()
}


# Get summaries and store these ##########################################

if (
  !file.exists(
    here::here("analysis/data/derived_data/beta-model-full-sum.RDS")
  )
) {
  bm1_sum <- seqwrap_summarise(
    readRDS(
      here::here("analysis/data/derived_data/beta-model-full.RDS")
    )
  )

  saveRDS(
    bm1_sum,
    here::here("analysis/data/derived_data/beta-model-full-sum.RDS")
  )
  rm(bm1_sum)
  gc()
}


if (
  !file.exists(
    here::here("analysis/data/derived_data/m-model-full-sum.RDS")
  )
) {
  mm1_sum <- seqwrap_summarise(
    readRDS(
      here::here("analysis/data/derived_data/m-model-full.RDS")
    )
  )

  saveRDS(
    mm1_sum,
    here::here("analysis/data/derived_data/m-model-full-sum.RDS")
  )
  rm(mm1_sum)
  gc()
}
