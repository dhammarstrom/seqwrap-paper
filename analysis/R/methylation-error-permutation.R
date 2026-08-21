# This script runs the permutation study on the methylation data set from
# seaborne et al..
# The resulting data is evaluated in case-study-methylation.qmd

# The script saves permutation results (one file per iteration)
# in analysis/data/derived_data/permutation.
# Analysis of the results are kept in the extended description of the case study.

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
# samping time-point. Evaluation requires thsi specific parameter from the model.
spike_time <- "3_loading"

# Check beta- and M-values are two representations of the same data
stopifnot(
  identical(rownames(m_vals), rownames(beta_vals)),
  identical(colnames(m_vals), colnames(beta_vals))
)

# We want to "test" the models on a subset of sites that represents low, mid
# and high beta values. row_id is assigned before any filtering so that it
# indexes rows of beta_vals / m_vals directly.
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

# plogis() saturates to exactly 0/1 near |M| = 37, which the beta model cannot
# accept, and the spike is what pushes extreme sites outwards. Check the worst
# case:
stopifnot(max(abs(m_vals)) + max(delta_grid) < 30)

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
# The four model arms are not all glmmTMB fits, and the accessors are not
# interchangeable. The eval_fun takes care of this, but returns different
# statistics per algorithm.
#
#   glmmTMB  -- S3 list. m$fit$convergence is 0 when the optimiser reports a
#               satisfactory solution; m$sdr$pdHess is TRUE when the Hessian is
#               positive definite (no singular fit, no unidentifiable
#               parameters). drop1(test = "Chisq") gives a chi-square LRT with
#               columns Df / AIC / LRT / Pr(>Chi).
#   merMod   -- S4 (lmerTest::lmer returns lmerModLmerTest).

# lme4::isSingular() is TRUE whenever a variance component sits on
# the zero boundary while glmmTMB reports only a non-positive-definite
# Hessian and flagged 0 of the
# same 120. Summing them into one "failure rate" would make the lmerTest arm
# look catastrophic when it is fitting fine. Report the two separately, and
# treat only `convergence != 0` and `is.na(pval)` as failures.
#
# `re_sd`, the estimated participant standard deviation, is recorded for every
# arm as the raw material for a comparable boundary diagnostic. Note that it is
# not directly threshold-able across arms either: on that same pilot the REML
# fits put 31 sites at exactly zero while the ML fits of the same data shrank
# to small-but-nonzero values (order 1e-6) at all but 2. Choose the threshold
# when analysing, and choose it per arm.
#
# Columns are therefore named for what they hold in general (`stat`, `test`)
# rather than for the chi-square case, and `ddf` is NA for the LRT arms.
# Columns are addressed BY NAME: positional indexing is what silently produced
# NAs when the merMod table turned out to have a different shape and no
# <none> row.
eval_fun <- function(m) {
  if (inherits(m, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(m)$cond$participant

    # A REML fit has NO valid likelihood-ratio test of a fixed effect: the
    # restricted likelihoods of the full and reduced models are taken over
    # different fixed-effect design matrices and are not comparable. drop1()
    # does not refuse -- it returns a NEGATIVE statistic and p = 1 (measured:
    # LRT = -10.53 for the beta arm, -7.51 for the Gaussian one). Those would
    # enter the error rates as guaranteed non-rejections and silently drag the
    # type I estimate towards zero, so the omnibus is left empty instead. The
    # REML arms contribute Wald contrasts only.
    reml <- isTRUE(m$modelInfo$REML)
    a <- if (reml) NULL else drop1(m, test = "Chisq")

    data.frame(
      convergence = m$fit$convergence,
      singular = !m$sdr$pdHess,
      re_sd = sqrt(as.numeric(vc[1, 1])),
      n = nrow(m$frame),
      test = if (reml) "none" else "LRT",
      df = if (reml) NA_real_ else a["time_permute", "Df"],
      ddf = NA_real_,
      stat = if (reml) NA_real_ else a["time_permute", "LRT"],
      pval = if (reml) NA_real_ else a["time_permute", "Pr(>Chi)"]
    )
  } else if (inherits(m, "merMod")) {
    a <- stats::anova(m, ddf = "Satterthwaite")
    cv <- m@optinfo$conv$opt
    vc <- lme4::VarCorr(m)$participant

    data.frame(
      # lme4 reports 0 for a satisfactory solution, as glmmTMB does.
      convergence = if (length(cv)) cv else NA_integer_,
      singular = lme4::isSingular(m),
      re_sd = sqrt(as.numeric(vc[1, 1])),
      n = nrow(m@frame),
      test = "Satt-F",
      df = a[["NumDF"]][1],
      ddf = a[["DenDF"]][1],
      stat = a[["F value"]][1],
      pval = a[["Pr(>F)"]][1]
    )
  } else {
    stop("unsupported model class: ", paste(class(m), collapse = "/"))
  }
}

# Wald summaries, one row per fixed effect.
#
# component = "cond" is passed explicitly rather than relying on it being the
# broom.mixed default for glmmTMB. Any arm carrying a dispersion or
# zero-inflation model has coefficients named time_permute* in that component
# too, and duplicated (model, term) pairs would silently corrupt the BH
# adjustment below -- that separation should not rest on a package setting.
#
# `df` is kept: for the lmerTest arm it holds the Satterthwaite denominator
# degrees of freedom, which is the whole point of that arm. glmmTMB fits have
# no such column, so it is filled with NA to keep the arms bind_rows-able.
summary_fun_wald <- function(m) {
  out <- broom.mixed::tidy(m, effects = "fixed", component = "cond")

  if (!"df" %in% names(out)) out$df <- NA_real_

  dplyr::select(out, effect, term, estimate, std.error, df, p.value)
}

summary_fun <- summary_fun_wald


#### Running the permutation loop #############################################

n_perm <- 500 #

# The results directory is VERSIONED, and the loop skips an iteration whose
# file already exists. analysis/data/derived_data/permutation holds the earlier
# run, which fitted two model arms only (beta, m) and stored the evaluation
# columns under the old names (pdHess, lrt). Writing the four-arm results into
# that directory would have skipped every iteration -- and had the files been
# removed instead, the two runs could not be told apart. Bump the suffix
# whenever the set of model arms or the eval_fun columns change.
out_dir <- here::here("analysis/data/derived_data/permutation_v2")
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
# lmerTest/lme4 fit the m_reml_satt arm, so they do work on the workers too.
# Satterthwaite degrees of freedom are computed by lmerTest itself and need no
# further package (Kenward-Roger would have pulled in pbkrtest).
.pkgs <- c(
  "seqwrap",
  "glmmTMB",
  "broom.mixed",
  "lmerTest",
  "lme4"
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
      # The design crosses SCALE (beta vs M) with ESTIMATOR (ML vs REML) and
      # REFERENCE DISTRIBUTION (asymptotic vs small-sample), giving seven
      # inferential cells from four fits:
      #
      #   arm            fit                          cells it supplies
      #   ------------   --------------------------   -----------------------
      #   beta           glmmTMB, beta, ML            beta ML Wald, beta ML LRT
      #   beta_reml      glmmTMB, beta, REML          beta REML Wald
      #   m              glmmTMB, gaussian, ML        m ML Wald, m ML LRT
      #   m_reml_satt    lmerTest::lmer (REML)        m REML Satterthwaite Wald
      #                                               + Satterthwaite F omnibus
      #
      # The seventh cell, "m REML Wald against a normal reference", needs no
      # fit of its own: glmmTMB(REML = TRUE) on the Gaussian model and
      # lmerTest::lmer return the same estimates and standard errors (checked
      # on this data: max absolute SE difference 2.9e-10), so it is recovered
      # downstream as 2 * pnorm(-|estimate / std.error|) from the
      # m_reml_satt summaries. Fitting it separately would cost a quarter more
      # compute for a numerically identical answer.
      #
      # Asymmetry worth noting when reading the results: the beta arms have no
      # small-sample omnibus available. REML removes the likelihood-ratio test
      # (see eval_fun), and there is no Satterthwaite or Kenward-Roger
      # machinery for a beta GLMM, so the only omnibus on the beta scale is
      # the asymptotic ML chi-square. That is a property of the available
      # tooling, not of this script.
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
        modelfun = lmerTest::lmer,
        data = m_subset,
        metadata = metapermute,
        samplename = "geo_accession",
        eval_fun = eval_fun,
        summary_fun = summary_fun,
        arguments = list(
          formula = y ~ time_permute + (1 | participant)
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
        tag(mm2_sum$summaries, "m_reml_satt", "summaries")
      ) |>
        # Keep the time contrasts
        dplyr::filter(base::grepl("^time", term)) |>
        # FDR is controlled ACROSS SITES within a term, not within a site
        dplyr::mutate(
          .by = c(model, term),
          qval = stats::p.adjust(p.value, method = "BH")
        ) |>
        dplyr::mutate(iter = i) |>
        dplyr::select(
          model,
          iter,
          target,
          term,
          estimate,
          std.error,
          # Satterthwaite denominator df for the m_reml_satt arm, NA elsewhere
          df,
          p.value,
          qval
        )

      eval_i <- dplyr::bind_rows(
        tag(bm1_sum$evaluations, "beta", "evaluations"),
        tag(bm2_sum$evaluations, "beta_reml", "evaluations"),
        tag(mm1_sum$evaluations, "m", "evaluations"),
        tag(mm2_sum$evaluations, "m_reml_satt", "evaluations")
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
