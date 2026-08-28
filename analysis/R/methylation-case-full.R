# This script fits the full data set from Seaborne et al. Raw data is
# downloaded in methylation-case-study-dataprep.R; here we proceed from the
# normalized gset.
#
# The model arms are the REML arms of the permutation study
# (analysis/R/methylation-error-permutation.R). Type I error rates from the
# permutation apply to the tests reported here:
#
# arm         fit                     tests reported
# ---------   ---------------------   ---------------------------
# m_reml      gaussian, REML          Wald-z, Wald-Satterthwaite
# beta_reml   beta_family, REML       Wald-z, Wald-Satterthwaite
#
# THE ML ARMS ARE NOT FITTED HERE. The permutation study is what compares the
# estimators, and it found the REML arms to be the calibrated ones; refitting
# the ML arms over every position would spend 7.2 hours reproducing an
# estimator this study has already rejected, and the shrinkage step
# (analysis/R/methylation-shrinkage.R) reads the REML arms only. The comparison
# lives in the permutation study, where it is answered at a hundredth of the
# cost.
#
# Two things differ from the permutation study by design:
#
# - No omnibus LRT. The permutation study establishes the type I error
# behaviour of the omnibus test; the full run reports per-contrast tests
# only, so eval_fun does NOT call drop1() and therefore refits nothing.
# drop1(test = "Chisq") refits the model once which is expensive.
#
# - Every arm here is read with Satterthwaite, because every arm here is REML.
# The permutation study is what fills the whole estimator x reference grid;
# ML+Satterthwaite was dropped from it on the grounds that it did not fix the
# type I error rate and cost a great deal, and plain ML is dropped from this
# script for the reason given above.
#
# Where Satterthwaite IS computed, both reference distributions come out of
# a single fit, in one pass of summary_fun (see summary_fun_wald below): the
# denominator df needs a second differentiation of the likelihood, which is
# not free, but it is paid once per fit, on the worker, and does not require
# the model to be fitted again.

## Packages and settings ###################################################
library(GEOquery) # For accessing the raw data
library(minfi) # For pre-processing of raw data
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(maxprobes) # For filtering cross-reactive probes
library(seqwrap)


# Detect number of cores. Every seqwrap() call below runs at this width, so
# the run adapts to whatever machine it is on. detectCores() is documented to
# return NA on systems it cannot interrogate, which would reach seqwrap() as
# `cores = NA` rather than as an error, so it is floored here.
cores <- parallel::detectCores()
if (is.na(cores) || cores < 1L) cores <- 1L


# The results directory is VERSIONED, in the same way and for the same reason
# as the permutation study: an arm is skipped when its file already exists, so
# writing a new set of arms, or a new set of summary/eval columns, into an
# existing directory would silently mix formats. Bump the suffix whenever the
# arms or the summary_fun / eval_fun columns change. Earlier runs:
#
#   (derived_data root)  two arms (beta, m) written as beta-model-full.RDS and
#                        m-model-full.RDS, default summary_fun (a single
#                        `p.value` column), eval_fun reporting pdHess. Those
#                        files are superseded by this directory.
#   full_v2              four glmmTMB arms (beta, beta_reml, m, m_reml), both
#                        Wald reference distributions per contrast, no LRT,
#                        fitted on the FUNCTIONALLY normalized arrays.
#   full_v3              the two REML arms only, same two Wald readings, on the
#                        QUANTILE normalized arrays. The normalization changed
#                        because the level-dependent shift funnorm leaves
#                        behind made the shrinkage step uninterpretable; see
#                        the note at `gset_file` in
#                        analysis/R/methylation-case-study-dataprep.R. Results
#                        from full_v2 and full_v3 are not comparable position
#                        by position and must not be pooled.
out_dir <- here::here("analysis/data/derived_data/full_v3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Quality control and normalization of data ##################################
#
# Moved to analysis/R/methylation-case-study-dataprep.R. The normalized gset is
# an input to both this script and the permutation study
# (analysis/R/methylation-error-permutation.R), so it is built once, upstream,
# rather than as a side effect of the script that fits the full models.

# Load normalized saved gset ################################################

gset_file <- here::here(
  "analysis/data/derived_data/seaborne-gset-quantile.RDS"
)

if (!file.exists(gset_file)) {
  stop(
    "seaborne-gset-quantile.RDS is missing. Run ",
    "analysis/R/methylation-case-study-dataprep.R first -- it downloads the ",
    "arrays if needed and writes the normalized, filtered gset.",
    call. = FALSE
  )
}

gset <- readRDS(gset_file)
metadata <- readRDS(
  here::here("analysis/data/derived_data/seaborne-metadata.RDS")
)

# Deriving M- and beta-values ##############################################
#
# M first, beta derived from it -- the reverse of what this script did under
# functional normalization, and the reversal is required rather than
# stylistic. preprocessQuantile() returns a GenomicRatioSet that stores
# M = log2(Meth / Unmeth) and no intensities at all, so the Illumina offset
# that used to hold beta off 0 and 1, getBeta(gset, offset = 100), has nothing
# left to act on: getBeta() derives beta from M by the inverse logit and
# ignores `offset` without complaint. Deriving beta from M explicitly is the
# same computation with the silence removed. See the normalization note in
# analysis/R/methylation-case-study-dataprep.R.
m_vals <- minfi::getM(gset)
beta_vals <- 2^m_vals / (1 + 2^m_vals)

# The inverse logit lands in the OPEN interval for every finite M, so no
# squeeze is needed. That is a property of the normalization rather than of
# this script, so it is checked here as well as where the gset is written --
# a beta of exactly 0 or 1 is a response the beta family cannot take, and it
# would surface hours into the run rather than now.
stopifnot(
  all(is.finite(m_vals)),
  all(beta_vals > 0 & beta_vals < 1)
)

# Clean up
# Remove gset as no longer used
rm(gset)
gc()


# Save data in suitable data frames ########################################
betadf <- data.frame(beta_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

mdf <- data.frame(m_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

# The two matrices hold the same data as the two data frames, and each is
# ~230 MB at this number of sites. The fits below run for hours with a results
# object of comparable size alive in memory, so the duplicates are removed here.
rm(beta_vals, m_vals)
gc()


# Evaluation function for the full setup ##################################
#
# One row per model fit: convergence diagnostics only. There is NO omnibus
# test here, and therefore no drop1() and no refitting, see above.
#
# m$fit$convergence is 0 when the optimiser reports a satisfactory solution.
# m$sdr$pdHess is TRUE when the Hessian is positive definite (no singular fit,
# no unidentifiable parameters); it is stored as `singular = !pdHess`, which
# is the spelling the permutation results use, so the two runs can be pooled.
# see e.g., https://stackoverflow.com/questions/79110546/glmmtmb-convergence-messages
# and https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html
#
# `re_sd`, the estimated participant standard deviation, is recorded for every
# fit; it is the quantity REML is chosen for, so it is what the ML/REML
# contrast is read on.
#
# `reml` is taken from the fit rather than from the arm name, so the arm
# labelling cannot silently disagree with what was actually fitted.
eval_fun <- function(m) {
  # Safety stop, the arms below are all glmmTMB fits
  if (!inherits(m, "glmmTMB")) {
    stop("unsupported model class: ", paste(class(m), collapse = "/"))
  }

  vc <- glmmTMB::VarCorr(m)$cond$participant

  data.frame(
    convergence = m$fit$convergence,
    singular = !m$sdr$pdHess,
    re_sd = sqrt(as.numeric(vc[1, 1])),
    n = nrow(m$frame),
    reml = isTRUE(m$modelInfo$REML)
  )
}


# Summary function ########################################################
#
# Kept deliberately IDENTICAL to summary_fun_wald in
# analysis/R/methylation-error-permutation.R: the permutation study measures
# the type I error rate of exactly these columns, so any divergence would make
# the error rates quoted for the full run rates of a different test.
#
# One row per fixed effect and two reference distributions per row:
#
# p_wald: the asymptotic Wald test, estimate / std.error against a normal
# p_satt: the same statistic against a t distribution with Satterthwaite
#           denominator degrees of freedom (`ddf`)
#
# Both come from coef(summary(.))$cond, i.e. from the SAME fit, so the only
# thing that differs between the two columns is the reference distribution.
# That is what makes the small-sample correction readable on its own,
# separately from the ML/REML contrast carried by the arms -- and it is why
# the two do not need two runs.
#
# summary(ddf = "satterthwaite") is a glmmTMB >= 1.1.13 feature and works for
# the beta family as well as the Gaussian one (confirmed in the permutation
# exp.), this allows the small-sample cell to be filled on both models
# (M and beta)
#
# The Satterthwaite computation needs a second differentiation of the
# likelihood and can fail on a degenerate fit where the plain Wald table is
# still returned. It is wrapped so that such a site loses `p_satt` only,
# rather than dropping out of every arm-by-test cell at once. Note that it can
# also SUCCEED and return a negative `ddf` on a singular fit, in which case
# p_satt comes back NaN -- so downstream aggregation must treat p_satt as
# possibly missing even where convergence == 0.
#
# `ddf` is stored per contrast, not just used, because it is also the
# diagnostic for this reading: on a site whose participant variance sits on
# the zero boundary the denominator df can collapse to well below the nominal
# residual df, and a p_satt from such a fit is conservative for a reason that
# has nothing to do with the small-sample correction working as intended.
summary_fun_wald <- function(m) {
  cf <- coef(summary(m))$cond

  out <- data.frame(
    term = rownames(cf),
    estimate = cf[, "Estimate"],
    std.error = cf[, "Std. Error"],
    p_wald = cf[, "Pr(>|z|)"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  sat <- tryCatch(
    coef(summary(m, ddf = "satterthwaite"))$cond,
    error = function(e) NULL
  )

  # Matched by term name, never by position.
  i <- if (is.null(sat)) NA_integer_ else match(out$term, rownames(sat))

  out$ddf <- if (is.null(sat)) NA_real_ else sat[i, "ddf"]
  out$p_satt <- if (is.null(sat)) NA_real_ else sat[i, "Pr(>|t|)"]

  out
}

# Model arms ################################################################
#
# SCALE (beta vs M) crossed with ESTIMATOR (ML vs REML). The response data
# frame differs between the scales; the estimator is the `REML` argument, and
# nothing else changes, so an arm is fully described by (data, arguments) plus
# the reference distributions it is read with (summary_fun).
#
# Arm names match the permutation study exactly.
arms <- list(
  list(
    name = "beta_reml",
    scale = "beta",
    data = betadf,
    summary_fun = summary_fun_wald,
    arguments = list(
      formula = y ~ time + (1 | participant),
      family = glmmTMB::beta_family(),
      REML = TRUE
    )
  ),
  list(
    name = "m_reml",
    scale = "M",
    data = mdf,
    summary_fun = summary_fun_wald,
    arguments = list(
      formula = y ~ time + (1 | participant),
      REML = TRUE
    )
  )
)

# Every arm here is REML and every arm is read with Satterthwaite, which is
# the whole point of dropping the ML arms. The arm table is where that could
# quietly drift -- a missing REML = TRUE would fit a different estimator under
# a name that claims otherwise, and the permutation study's error rates would
# then be quoted for a test that was never run. Assert it rather than trusting
# the list above.
stopifnot(
  vapply(arms, function(a) isTRUE(a$arguments$REML), logical(1)),
  vapply(
    arms,
    function(a) identical(a$summary_fun, summary_fun_wald),
    logical(1)
  ),
  vapply(arms, function(a) grepl("_reml$", a$name), logical(1))
)


# Fit the arms ##############################################################
#
# One seqwrap() call per arm, each timed and cached on its own. A completed arm
# is skipped, so an interrupted run resumes at the arm it stopped in rather
# than refitting the ones already done.
#
# Only the SUMMARISED results are written. seqwrap() is called with
# return_models = FALSE, so its result object carries the same summaries that
# seqwrap_summarise() returns, only in unbound per-target form -- saving both
# doubled the space on disk for no extra information. Everything the analysis
# reads is carried across explicitly below: summaries, evaluations, errors,
# the elapsed time and the target counts.
#
# rm()/gc() between arms matters: a results object at this number of sites is
# on the order of 10^2 MB, and four of them alive at once is what would push
# this run into swap.
for (arm in arms) {
  arm_file <- file.path(out_dir, sprintf("%s.RDS", arm$name))
  meta_file <- file.path(out_dir, sprintf("%s-meta.RDS", arm$name))

  if (file.exists(arm_file)) {
    message("skipping arm ", arm$name, " -- ", basename(arm_file), " exists")
    next
  }

  message("fitting arm ", arm$name, " on ", cores, " cores")

  container <- seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    data = arm$data,
    metadata = metadata,
    samplename = "geo_accession",
    eval_fun = eval_fun,
    summary_fun = arm$summary_fun,
    arguments = arm$arguments
  )

  results <- seqwrap(
    container,
    return_models = FALSE,
    cores = cores,
    #   subset = 1:10,
    verbose = FALSE
  )

  summarised <- seqwrap_summarise(results, verbose = FALSE)

  # seqwrap_summarise() OMITS an element entirely when every target failed at
  # that stage -- it does not return an empty data frame. Catching it here
  # names the arm, rather than letting a NULL surface hours later in the
  # analysis.
  if (is.null(summarised$summaries) || !NROW(summarised$summaries)) {
    stop(
      "arm ",
      arm$name,
      " returned no summaries. Every fit failed at that stage -- check that ",
      "eval_fun/summary_fun support this model class.",
      call. = FALSE
    )
  }

  # The timing is the reason each arm is its own seqwrap() call: elapsed_time
  # is measured across a single call, so one call per arm is what makes the
  # per-arm cost readable. It is stored both in the arm file and in a small
  # sidecar, so the timings can be read back without loading ~10^2 MB per arm.
  meta <- data.frame(
    model = arm$name,
    scale = arm$scale,
    reml = isTRUE(arm$arguments$REML),
    # Which reference distributions this arm actually paid for. Recorded so
    # that the timings stay interpretable: the Satterthwaite arms are not
    # slower because REML is slower -- REML fits are the FASTER ones here.
    satterthwaite = identical(arm$summary_fun, summary_fun_wald),
    cores = cores,
    k = results@k,
    n = results@n,
    n_summarised = length(unique(summarised$summaries$target)),
    user = results@elapsed_time[["user.self"]],
    system = results@elapsed_time[["sys.self"]],
    elapsed = results@elapsed_time[["elapsed"]],
    n_conditions = NROW(summarised$errors),
    stringsAsFactors = FALSE
  )

  saveRDS(
    list(
      model = arm$name,
      arguments = results@call_arguments,
      engine = results@call_engine,
      elapsed_time = results@elapsed_time,
      k = results@k,
      n = results@n,
      summaries = summarised$summaries,
      evaluations = summarised$evaluations,
      errors = summarised$errors,
      meta = meta
    ),
    arm_file
  )

  saveRDS(meta, meta_file)

  message(
    "arm ",
    arm$name,
    " done in ",
    round(meta$elapsed / 60 / 60, 2),
    " hours"
  )

  # Clean-up ##
  rm(container, results, summarised, meta)
  gc()
}


# Collected timings #########################################################
#
# This is the table the computational-time section of the case
# study.
meta_files <- file.path(
  out_dir,
  sprintf("%s-meta.RDS", vapply(arms, `[[`, character(1), "name"))
)

if (all(file.exists(meta_files))) {
  timing <- dplyr::bind_rows(lapply(meta_files, readRDS)) |>
    dplyr::mutate(
      estimator = dplyr::if_else(reml, "REML", "ML"),
      hours = elapsed / 60 / 60
    )

  saveRDS(timing, file.path(out_dir, "timing.RDS"))

  print(timing[, c("model", "scale", "estimator", "k", "n", "hours")])
}
