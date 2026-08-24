# This script runs the permutation study on the methylation data set from
# seaborne et al. The resulting data is evaluated in case-study-methylation.qmd
# and in the main paper source.

# The script saves permutation results (one file per iteration) in the
# versioned results directory set by `out_dir` below
# (analysis/data/derived_data/permutation_v3).
# Analysis of the results are kept in the extended description
# of the case study.

# Load packages ###########################
library(seqwrap)
library(tidyverse)
library(minfi)


# Load data ##############################
# Normalized gset
gset <- readRDS(
  here::here(
    "analysis/data/derived_data/seaborne-gset-normalized.RDS"
  )
)
# Metadata
metadata <- readRDS(
  here::here("analysis/data/derived_data/seaborne-metadata.RDS")
)


# Calculate the beta-values with the Illumina-style offset (c = 100)
beta_vals <- getBeta(gset, offset = 100)
# Derive M-values from the calculated beta values.
m_vals <- log2(beta_vals / (1 - beta_vals))


# Simulation settings #####################################

# Over all seed for design etc.
set.seed(1)

# Set the number of cores to use - default to detect cores
CORES <- parallel::detectCores()

# Sites sampled per beta-value stratum
size <- 200

# Effect sizes (in M-value units) multiplied with a sign, and the number of
# spiked sites allocated to each (delta x sign) cell within every stratum.
delta_grid <- c(0.1, 0.2, 0.4, 0.8)
cells <- expand.grid(delta = delta_grid, sign = c(1, -1))
n_per_cell <- 5

# Sites spiked per stratum, and the resulting proportion of true alternatives
n_spike <- nrow(cells) * n_per_cell
pi1 <- n_spike / size


# delta_grid contains the effect magnitudes and the direction (sign) of the
# effect. See below for assign_spike().

# The delta grid must...
stopifnot(
  # not contain zero entries
  all(delta_grid > 0),
  # not contain duplicate
  !anyDuplicated(delta_grid),
  # Total size should be > than the number spikes
  size > n_spike
)

# The time-point that receives the simulated effect is the 3_loading
# samping time-point. Evaluation requires thsi specific parameter
# from the model.
spike_time <- "3_loading"

# Check beta- and M-values are two representations of the same data
stopifnot(
  identical(rownames(m_vals), rownames(beta_vals)),
  identical(colnames(m_vals), colnames(beta_vals))
)

# We want to "test" the models on a subset of sites that represents low, mid
# and high beta values. row_id is assigned before any filtering so that it
# indexes rows of beta_vals / m_vals directly.
# The subset is fixed in permutations - The seed above determines
# the sites selected here.

subset <- data.frame(
  id = row.names(beta_vals),
  betacat = cut(
    rowMeans(beta_vals),
    breaks = seq(0, 1, by = 0.2),
    labels = c("xlow", "low", "mid", "high", "xhigh"),
    include.lowest = TRUE
  )
) |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  dplyr::filter(!is.na(betacat)) |>
  dplyr::group_by(betacat) |>
  dplyr::slice_sample(n = size) |>
  dplyr::ungroup()

# Site-by-sample M-value matrix for the sampled sites, in `subset` row order.
# Every positional operation below relies on this alignment.
m_sub <- m_vals[subset$row_id, , drop = FALSE]
stopifnot(identical(rownames(m_sub), subset$id))

# plogis() saturates (in practice) to exactly 0/1 near |M| = 17, which
# breaks the beta model. The added effect adds to extreme sites.
# Here we check the worst case:
stopifnot(max(abs(m_vals)) + max(delta_grid) < 17)


# Permute time labels within participant. This keeps each participants mean
# and the within-participant correlation, and removes the association between
# the response and the time-point.
permute_time <- function(md) {
  md |>
    dplyr::group_by(participant) |>
    dplyr::mutate(time_permute = sample(time)) |>
    dplyr::ungroup() |>
    dplyr::select(geo_accession, participant, time, time_permute)
}

# Allocate spiked sites to (delta x sign) cells, balanced within stratum.
# Sites that are not spiked get delta = 0 and are the true nulls.
#
# `delta` is returned SIGNED, i.e. it is the effect actually added to the
# M-values, so it can be used directly downstream. `delta_abs` keeps the
# magnitude for summaries that pool the two directions.

# sub, the subset data frame,
# cells, cells in the effect size grid
# n_per_cell, number of strata (5) per effect size grid cell.
assign_spike <- function(sub, cells, n_per_cell) {
  # One code per spiked site: cell 1 repeated n_per_cell times, then cell 2,
  # and so on. Codes are 1..nrow(cells); 0 marks an unspiked site.
  codes <- rep(seq_len(nrow(cells)), each = n_per_cell)

  # Signed effect per cell, with a leading 0 so that code 0 -> no effect
  cell_delta <- c(0, cells$delta * cells$sign)

  sub |>
    dplyr::group_by(betacat) |>
    # Pad the codes with zeros up to the size of the stratum and shuffle:
    # every cell gets exactly n_per_cell sites, allocated at random.
    dplyr::mutate(
      cell_id = sample(
        c(codes, rep(0L, dplyr::n() - length(codes)))
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(delta = cell_delta[cell_id + 1L], delta_abs = abs(delta))
}


# Per-fit evaluation: convergence diagnostics plus an omnibus test of the time
# factor. The omnibus test has one row per model fit and it is therefore
# placed in eval_fun rather than in summary_fun.
#
# All four arms are fitted with glmmTMB, so the diagnostics are directly
# comparable across arms. The earlier lmerTest arm is gone: glmmTMB computes
# Satterthwaite degrees of freedom itself (see summary_fun below), which also
# makes the small-sample reference distribution available on the beta scale,
# where lmerTest could not follow.
#
# m$fit$convergence is 0 when the optimiser reports a satisfactory solution;
# m$sdr$pdHess is TRUE when the Hessian is positive definite (no singular fit,
# no unidentifiable parameters). The two are reported separately, and only
# `convergence != 0` and an unusable p-value count as failures.
#
# For the omnibus test, drop1(test = "Chisq") gives a chi-square LRT with
# columns Df / AIC / LRT / Pr(>Chi).
#
# `re_sd`, the estimated participant standard deviation, is recorded for every
# model for potential comparisons between models.
#
# Columns are named for what they hold in general (`stat`, `test`) rather than
# for the chi-square case. Columns are addressed BY NAME: positional indexing
# is what silently produced NAs when a table turned out to have a different
# shape and no <none> row.
eval_fun <- function(m) {
  # Safety stop, the arms below are all glmmTMB fits
  if (!inherits(m, "glmmTMB")) {
    stop("unsupported model class: ", paste(class(m), collapse = "/"))
  }

  vc <- glmmTMB::VarCorr(m)$cond$participant

  # A REML fit has no valid likelihood-ratio test of a fixed effect: the
  # restricted likelihoods of two different fixed-effect designs are not
  # comparable. The omnibus is left empty for the REML arms, which therefore
  # contribute Wald contrasts only.

  # Check if the fit is REML
  reml <- isTRUE(m$modelInfo$REML)
  # Save the drop1 output
  a <- if (reml) NULL else drop1(m, test = "Chisq")

  data.frame(
    convergence = m$fit$convergence,
    singular = !m$sdr$pdHess,
    re_sd = sqrt(as.numeric(vc[1, 1])),
    n = nrow(m$frame),
    reml = reml,
    test = if (reml) "none" else "LRT",
    df = if (reml) NA_real_ else a["time_permute", "Df"],
    stat = if (reml) NA_real_ else a["time_permute", "LRT"],
    pval = if (reml) NA_real_ else a["time_permute", "Pr(>Chi)"]
  )
}

# Wald summaries, one row per fixed effect and two reference distributions per
# row:
#
#   p_wald  the asymptotic Wald test, estimate / std.error against a normal
#   p_satt  the same statistic against a t distribution with Satterthwaite
#           denominator degrees of freedom (`ddf`)
#
# Both come from coef(summary(.))$cond, i.e. from the same fit, so the ONLY
# thing that differs between the two columns is the reference distribution.
# That is what makes the small-sample correction readable on its own,
# separately from the ML/REML contrast carried by the arms.
#
# summary(ddf = "satterthwaite") is a glmmTMB >= 1.1.11 feature and works for
# the beta family as well as the Gaussian one, which is what allows the
# small-sample cell to be filled on both scales.
#
# The `cond` component is taken explicitly rather than relying on a default:
# any arm carrying a dispersion or zero-inflation model has coefficients named
# time_permute* in that component too, and duplicated (model, term) pairs
# would silently corrupt the BH adjustment applied downstream.
#
# The Satterthwaite computation needs a second differentiation of the
# likelihood and can fail on a degenerate fit where the plain Wald table is
# still returned. It is wrapped so that such a site loses `p_satt` only,
# rather than dropping out of every arm-by-test cell at once.
#
# `ddf` is stored per contrast, not just used, because it is also the
# diagnostic for this arm: on a site whose participant variance sits on the
# zero boundary the denominator df can collapse to well below the nominal
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

summary_fun <- summary_fun_wald


#### Running the permutation loop #############################################

n_perm <- 500 #

# The results directory is VERSIONED, and the loop skips an iteration whose
# file already exists. Writing a new set of arms into an existing directory
# would skip every iteration -- and had the files been removed instead, the
# runs could not be told apart. Bump the suffix whenever the set of model arms
# or the eval_fun / summary_fun columns change. Earlier runs:
#
#   permutation      two arms (beta, m), old eval column names (pdHess, lrt)
#   permutation_v2   four arms (beta, beta_reml, m, m_reml_satt), where the
#                    Satterthwaite arm was an lmerTest fit on the M scale only
#   permutation_v3   four glmmTMB arms (beta, beta_reml, m, m_reml) with both
#                    Wald reference distributions stored per contrast
out_dir <- here::here("analysis/data/derived_data/permutation_v3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Worker start-up cost.
#
# seqwrap() creates a PSOCK cluster and destroys it again on every call, and
# each worker runs the project .Rprofile -- so renv activates once per worker,
# per call. Measured on this machine that is ~19 s of fixed cost per seqwrap()
# call at cores = 12, against ~22 s of actual fitting for 1000 sites; the loop
# makes four calls per iteration, so most of the run would be start-up.
#
# Pointing R_PROFILE_USER at an empty file makes the workers skip .Rprofile,
# and pushing .libPaths() through R_LIBS keeps them on the same renv library so
# package versions are unchanged. This takes the fixed cost to ~4 s per call.
.env_saved <- Sys.getenv(c("R_PROFILE_USER", "R_LIBS"), unset = NA)
.empty_profile <- tempfile(fileext = ".R")
file.create(.empty_profile)
Sys.setenv(
  R_PROFILE_USER = .empty_profile,
  R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)
)

# pblapply prints a progress bar on every seqwrap() call; noise in a loop.
.pbo <- pbapply::pboptions(type = "none")

# Because the workers no longer run .Rprofile, renv does not hide the base R
# library from them: they get the renv library and sandbox first (identical to
# this session) with the system library appended as a fallback. Packages that
# renv manages therefore resolve identically -- but one that is MISSING from the
# renv library would fail here and be silently taken from the system library on
# a worker. Assert that the packages doing the work resolve to the same place.
#
# normalizePath() matters: for a package that is *attached* in this session
# find.package() reports the path it was loaded from, and renv's library entries
# are junctions into the renv cache -- so the parent would say ".../renv/cache/
# .../seqwrap" where a worker, resolving through .libPaths(), says ".../renv/
# library/.../seqwrap". Same installation, two spellings. Comparing the resolved
# directory removes that false alarm while still catching a genuine fallback to
# the system library.
# All four arms are glmmTMB fits, and the Satterthwaite degrees of freedom are
# computed by glmmTMB itself (>= 1.1.11), so no further package is needed on
# the workers -- Kenward-Roger would have pulled in pbkrtest, and the earlier
# lmerTest arm is gone.
.pkgs <- c(
  "seqwrap",
  "glmmTMB"
)
.stamp <- function(pk)
  vapply(
    pk,
    function(p)
      paste(
        normalizePath(find.package(p), winslash = "/", mustWork = FALSE),
        utils::packageVersion(p)
      ),
    character(1)
  )
.parent <- .stamp(.pkgs)
.cl <- parallel::makeCluster(2)
.worker <- parallel::clusterCall(.cl, .stamp, .pkgs)
parallel::stopCluster(.cl)
stopifnot(all(vapply(.worker, identical, logical(1), .parent)))

# tryCatch(finally =) rather than local(): the restore below runs even if the
# loop aborts, while assignments inside the loop stay visible in the global
# environment for inspection after a failure.
tryCatch(
  {
    for (i in seq_len(n_perm)) {
      iter_file <- file.path(out_dir, sprintf("perm_%04d.RDS", i))
      # IF the file is already there, go to next iteration
      if (file.exists(iter_file)) next

      set.seed(1000 + i)

      ## (1) Permute the time labels within participant
      metapermute <- permute_time(metadata)
      stopifnot(setequal(colnames(m_sub), metapermute$geo_accession))

      ## (2) Draw a spike assignment
      spike_design <- assign_spike(subset, cells, n_per_cell)
      stopifnot(identical(spike_design$id, subset$id))

      ## (3) Samples that receive the simulated effect
      # Matched by name, never by position: the column order of m_sub and the
      # row order of metapermute are not guaranteed to agree.
      x <- as.integer(
        metapermute$time_permute[
          match(colnames(m_sub), metapermute$geo_accession)
        ] ==
          spike_time
      )
      stopifnot(!anyNA(x), sum(x) > 0)

      ## (4) Add the effect on the M scale, derive beta from it ------------
      # delta is already signed, and is zero for unspiked rows, while x is zero
      # for unspiked columns -- so outer() places the shift on the intersection
      # only.
      s <- spike_design$delta

      m_spiked <- m_sub + outer(s, x)
      beta_spiked <- stats::plogis(m_spiked * log(2))

      # plogis saturates to exactly 0/1 around |M| = 37, which the beta model
      # cannot accept. The spike is what pushes extreme sites outwards.
      stopifnot(max(abs(m_spiked)) < 30, all(beta_spiked > 0 & beta_spiked < 1))

      m_subset <- data.frame(id = subset$id, m_spiked, check.names = FALSE)
      beta_subset <- data.frame(
        id = subset$id,
        beta_spiked,
        check.names = FALSE
      )

      ## (5) Fit the model arms to the same permuted, spiked data ---------
      #
      # The design crosses SCALE (beta vs M) with ESTIMATOR (ML vs REML), and
      # each fit is then read with every REFERENCE DISTRIBUTION available to
      # it. Four fits, all glmmTMB, and ten inferential cells:
      #
      #   arm         fit                     cells it supplies
      #   ---------   ---------------------   ---------------------------
      #   m           gaussian, ML            Wald-z, Wald-Satterthwaite,
      #                                       LRT (omnibus)
      #   m_reml      gaussian, REML          Wald-z, Wald-Satterthwaite
      #   beta        beta_family, ML         Wald-z, Wald-Satterthwaite,
      #                                       LRT (omnibus)
      #   beta_reml   beta_family, REML       Wald-z, Wald-Satterthwaite
      #
      # The two Wald cells per arm are two readings of ONE fit and come out of
      # summary_fun together (p_wald, p_satt), so the reference distribution
      # can be varied with the estimator held fixed, and vice versa. The LRT
      # is an omnibus test of the whole time factor and has one row per fit,
      # so it lives in eval_fun; REML supplies no LRT (see eval_fun).
      #
      # This replaces the earlier lmerTest arm. glmmTMB (>= 1.1.11) computes
      # Satterthwaite degrees of freedom itself, for the beta family as well
      # as the Gaussian one, so the small-sample correction is now available
      # on both scales rather than on the M scale alone -- and every cell
      # comes from the same fitting machinery, which removes lme4-vs-glmmTMB
      # differences (convergence codes, singularity flags) from the
      # comparison. On the Gaussian model the two implementations agree
      # anyway: glmmTMB(REML = TRUE) and lmerTest::lmer returned the same
      # estimates and standard errors on this data, to 2.9e-10.
      #
      # A fifth arm modelling dispersion by time (dispformula = ~time_permute)
      # was tried and dropped. Two measurements against this same permutation
      # null ruled it out. First, there is nothing to model: the LRT of
      # dispformula ~time against ~1 rejects at 23.9% on the real labels but
      # 33.4% on permuted labels, so the real data shows LESS apparent
      # heteroscedasticity than a null with none by construction. Second, it
      # makes the problem this study is about worse -- paired on the same
      # sites, Wald type I at alpha = 0.05 rose from 7.8% to 10.6%, because
      # four extra precision parameters estimated from 7-8 observations each
      # give noisy standard errors that a normal reference does not correct
      # for. The idea is sound in principle (it is the Welch situation), but
      # the small-sample reference distribution it needs does not exist for a
      # beta GLMM.
      bm1 <- seqwrap_compose(
        modelfun = glmmTMB::glmmTMB,
        data = beta_subset,
        metadata = metapermute,
        samplename = "geo_accession",
        eval_fun = eval_fun,
        summary_fun = summary_fun,
        arguments = list(
          formula = y ~ time_permute + (1 | participant),
          family = glmmTMB::beta_family()
        )
      )

      bm2 <- seqwrap_compose(
        modelfun = glmmTMB::glmmTMB,
        data = beta_subset,
        metadata = metapermute,
        samplename = "geo_accession",
        eval_fun = eval_fun,
        summary_fun = summary_fun,
        arguments = list(
          formula = y ~ time_permute + (1 | participant),
          family = glmmTMB::beta_family(),
          REML = TRUE
        )
      )

      mm1 <- seqwrap_compose(
        modelfun = glmmTMB::glmmTMB,
        data = m_subset,
        metadata = metapermute,
        samplename = "geo_accession",
        eval_fun = eval_fun,
        summary_fun = summary_fun,
        arguments = list(
          formula = y ~ time_permute + (1 | participant)
        )
      )

      mm2 <- seqwrap_compose(
        modelfun = glmmTMB::glmmTMB,
        data = m_subset,
        metadata = metapermute,
        samplename = "geo_accession",
        eval_fun = eval_fun,
        summary_fun = summary_fun,
        arguments = list(
          formula = y ~ time_permute + (1 | participant),
          REML = TRUE
        )
      )

      bm1_results <- seqwrap(
        bm1,
        return_models = FALSE,
        cores = CORES,
        verbose = FALSE
      )

      bm2_results <- seqwrap(
        bm2,
        return_models = FALSE,
        cores = CORES,
        verbose = FALSE
      )

      mm1_results <- seqwrap(
        mm1,
        return_models = FALSE,
        cores = CORES,
        verbose = FALSE
      )

      mm2_results <- seqwrap(
        mm2,
        return_models = FALSE,
        cores = CORES,
        verbose = FALSE
      )

      bm1_sum <- seqwrap_summarise(bm1_results, verbose = FALSE)
      bm2_sum <- seqwrap_summarise(bm2_results, verbose = FALSE)
      mm1_sum <- seqwrap_summarise(mm1_results, verbose = FALSE)
      mm2_sum <- seqwrap_summarise(mm2_results, verbose = FALSE)

      ## (6) Collect -------------------------------------------------------
      #
      # seqwrap_summarise() OMITS an element entirely when every target failed
      # at that stage -- it does not return an empty data frame. A model arm
      # whose eval_fun raises on every fit therefore arrives here as NULL, and
      # `NULL |> mutate()` fails with "no applicable method for 'mutate'
      # applied to an object of class NULL", several hundred model fits after
      # the actual cause. tag() turns that into a message naming the arm and
      # the stage.
      tag <- function(x, arm, slot) {
        if (is.null(x) || !NROW(x)) {
          stop(
            "model arm '",
            arm,
            "' returned no ",
            slot,
            ". ",
            "Every fit failed at that stage -- check that eval_fun/summary_fun ",
            "support this model class.",
            call. = FALSE
          )
        }
        dplyr::mutate(x, model = arm)
      }

      res_i <- dplyr::bind_rows(
        tag(bm1_sum$summaries, "beta", "summaries"),
        tag(bm2_sum$summaries, "beta_reml", "summaries"),
        tag(mm1_sum$summaries, "m", "summaries"),
        tag(mm2_sum$summaries, "m_reml", "summaries")
      ) |>
        # Keep the time contrasts
        dplyr::filter(base::grepl("^time", term)) |>
        # FDR is controlled ACROSS SITES within a term, not within a site.
        # Both reference distributions are adjusted, and separately: they are
        # different p-values, so the BH ranking within a (model, term) is not
        # the same in general and one q-value cannot stand for the other.
        dplyr::mutate(
          .by = c(model, term),
          qval_wald = stats::p.adjust(p_wald, method = "BH"),
          qval_satt = stats::p.adjust(p_satt, method = "BH")
        ) |>
        dplyr::mutate(iter = i) |>
        dplyr::select(
          model,
          iter,
          target,
          term,
          estimate,
          std.error,
          # Asymptotic Wald, normal reference
          p_wald,
          qval_wald,
          # Wald against t with Satterthwaite denominator df
          ddf,
          p_satt,
          qval_satt
        )

      eval_i <- dplyr::bind_rows(
        tag(bm1_sum$evaluations, "beta", "evaluations"),
        tag(bm2_sum$evaluations, "beta_reml", "evaluations"),
        tag(mm1_sum$evaluations, "m", "evaluations"),
        tag(mm2_sum$evaluations, "m_reml", "evaluations")
      ) |>
        dplyr::mutate(iter = i)

      saveRDS(
        list(
          iter = i,
          summaries = res_i,
          evaluations = eval_i,
          design = spike_design |>
            dplyr::select(id, betacat, delta, delta_abs) |>
            dplyr::mutate(iter = i),
          metapermute = metapermute
        ),
        iter_file
      )
    }
  },
  finally = {
    # Runs whether the loop completed or aborted, so the modified environment
    # never leaks into the rest of the document.
    pbapply::pboptions(.pbo)
    for (v in names(.env_saved)) {
      if (is.na(.env_saved[[v]])) Sys.unsetenv(v) else
        do.call(Sys.setenv, stats::setNames(list(.env_saved[[v]]), v))
    }
    unlink(.empty_profile)
  }
)
