# This script runs the permutation study on the methylation data set
# seaborne et al..
# The resulting data is evaluated in case-study-methylation.qmd

# The script saves permutation results (one file per iteration)
# in analysis/data/derived_data/permutation.
# Analysis of the results are kept in the extended description of the case study.

# Load packages

library(seqwrap)
library(tidyverse)
library(minfi)


# Load data
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

# beta- and M-values are two representations of the same data
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


# Convergence information for each fit.
# m$fit$convergence is 0 when the optimiser reports a satisfactory solution.
# m$sdr$pdHess is TRUE when the Hessian is positive definite (no singular fit,
# no unidentifiable parameters). We want 0 and TRUE.
eval_fun <- function(m) {
  data.frame(
    convergence = m$fit$convergence,
    pdHess = m$sdr$pdHess,
    n = nrow(m$frame)
  )
}

# Wald summaries only: one model fit per site. We will defualt to this, but
# see below for an omnibus test.

summary_fun_wald <- function(m) {
  broom.mixed::tidy(m) |>
    dplyr::select(effect, term, estimate, std.error, p.value)
}

# Wald summaries plus an omnibus likelihood-ratio test for the time factor.
# This refits a reduced model at every site and therefore doubles the number of
# model fits. The reduced model is built from the fitted model frame rather than
# with update(), so that it does not depend on objects in the calling
# environment -- update() is not reliable inside socket-cluster workers.
summary_fun_lrt <- function(m) {
  out <- broom.mixed::tidy(m) |>
    dplyr::select(effect, term, estimate, std.error, p.value)

  p_lrt <- tryCatch(
    {
      m0 <- glmmTMB::glmmTMB(
        y ~ (1 | participant),
        data = m$frame,
        family = stats::family(m)
      )
      stats::anova(m0, m)[["Pr(>Chisq)"]][2]
    },
    error = function(e) NA_real_
  )

  dplyr::bind_rows(
    out,
    data.frame(
      effect = "omnibus",
      term = "time_omnibus",
      estimate = NA_real_,
      std.error = NA_real_,
      p.value = p_lrt
    )
  )
}

# TRUE doubles the number of model fits. Keep FALSE for the pilot.
use_lrt <- FALSE
summary_fun <- if (use_lrt) summary_fun_lrt else summary_fun_wald


#### Running the permutation loop #############################################

n_perm <- 250 #
out_dir <- here::here("analysis/data/derived_data/permutation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Worker start-up cost.
#
# seqwrap() creates a PSOCK cluster and destroys it again on every call, and
# each worker runs the project .Rprofile -- so renv activates once per worker,
# per call. Measured on this machine that is ~19 s of fixed cost per seqwrap()
# call at cores = 12, against ~22 s of actual fitting for 1000 sites; the loop
# makes two calls per iteration, so most of the run would be start-up.
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
.pkgs <- c("seqwrap", "glmmTMB", "broom.mixed")
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

      ## (5) Fit both models to the same permuted, spiked data -------------
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

      bm1_results <- seqwrap(
        bm1,
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

      bm1_sum <- seqwrap_summarise(bm1_results, verbose = FALSE)
      mm1_sum <- seqwrap_summarise(mm1_results, verbose = FALSE)

      ## (6) Collect -------------------------------------------------------
      res_i <- dplyr::bind_rows(
        bm1_sum$summaries |> dplyr::mutate(model = "beta"),
        mm1_sum$summaries |> dplyr::mutate(model = "m")
      ) |>
        # Keep the time contrasts (and the omnibus test when enabled)
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
          p.value,
          qval
        )

      eval_i <- dplyr::bind_rows(
        bm1_sum$evaluations |> dplyr::mutate(model = "beta"),
        mm1_sum$evaluations |> dplyr::mutate(model = "m")
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
