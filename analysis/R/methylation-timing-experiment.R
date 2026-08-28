# Timing experiment: what does seqwrap() cost, and what does the cost scale
# with?
#
# The full run (analysis/R/methylation-case-full.R) fits 770 441 CpG positions
# per arm and takes between ~3 and ~10.6 hours per arm on 16 cores. That single
# number per arm says what the run cost, but not what drives the cost. This
# script measures the same four arms on sub-samples of the SAME data, so the
# cost can be read as a function of the three things that are varied here:
#
#   block "size"        k, the number of targets fitted, at fixed model and
#                       fixed sample size. The claim to be tested is that
#                       elapsed time is linear in k -- i.e. seqwrap adds a
#                       fixed overhead per call and a constant cost per target,
#                       and nothing that grows faster than k.
#
#   block "complexity"  the model, at fixed k and fixed sample size. A ladder
#                       of formulas from a fixed-effects-only model to two
#                       random-effect terms, so per-target cost can be read
#                       against the number of estimated parameters.
#
#   block "samples"     n, the number of arrays per target, at fixed k and
#                       fixed model. Participants are dropped whole, so the
#                       repeated-measures structure is kept.
#
#   block "anchor"      one fixed cell (see ANCHOR below) re-run at intervals
#                       throughout the experiment. It measures nothing about
#                       seqwrap; it measures the MACHINE -- thermal throttling,
#                       background load, anything that drifts over the ~2 hours
#                       the experiment takes. A trend in the anchor invalidates
#                       comparisons between cells run early and cells run late,
#                       and the run order is randomised (see below) so that
#                       such a drift shows up as noise rather than as a slope.
#
# What makes the extrapolation legitimate. The timed region is one seqwrap()
# call, and the work inside it is the work the full run did: same data, same
# `arguments`, same `summary_fun` (Satterthwaite on the REML arms is a real
# cost and is paid here too), same `eval_fun`, same `return_models = FALSE`,
# same `cores`, and the default chunk size, which is itself a function of k and
# the core count (ceiling(k / (4 * cores)), capped at 2000) and so keeps the
# parallel granularity roughly constant as k grows. The one difference is that
# the targets are a random sub-sample of the sites, which is what the size axis
# requires. The fitted intercept and slope are therefore in the same units as
# the full run, and the check at the end of the script is exactly that: predict
# the full run from the small runs and compare with the observed hours in
# full_v2/timing.RDS.
#
# Two design choices that are there to keep the size axis clean:
#
# - Nested sub-samples. Within a replicate, the k = 100 targets are the first
#   100 of the k = 4000 draw, the k = 250 targets the first 250, and so on. The
#   sites therefore do not change composition between the k levels of a
#   replicate: a site that is expensive to fit (slow convergence, boundary
#   variance) is present at every k above the point it enters. Composition
#   still varies BETWEEN replicates, which is what the replicates are for.
#
# - The data are sliced BEFORE seqwrap_compose(), and seqwrap()'s own `subset`
#   argument is not used. `subset` is applied inside the timed region, so
#   subsetting a 770 441-row data frame there would add a fixed cost of the
#   order of the per-call overhead this experiment is trying to measure -- and
#   the full run, which passed no subset, never paid it.
#
# Running order is randomised over all cells (one fixed seed), so that any
# drift in machine speed is not aliased onto the design. The machine should be
# otherwise idle: this measures wall-clock time, and a browser is a covariate.
#
# Runtime: roughly 1.5-2 hours with the settings below, on the same 16 cores
# the full run used. Set QUICK <- TRUE for a few-minute smoke test of the
# harness that produces a design of the same shape at k <= 250.
#
# Results are written one file per cell to the versioned directory set by
# `out_dir` below; a completed cell is skipped, so an interrupted run resumes.


## Packages and settings ###################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(minfi)
library(seqwrap)

# A smoke test of the harness rather than the experiment: same blocks, same
# code path, k small enough to finish in minutes. Timings from a QUICK run are
# NOT the experiment -- at these k the per-call overhead dominates.
QUICK <- FALSE

# Cores. Kept identical to the full run, since the quantity being extrapolated
# is the cost AT THAT CORE COUNT: the per-target slope is a function of how
# many workers share the targets, and the intercept is largely cluster
# start-up, which grows with the number of workers.
CORES <- parallel::detectCores()

# Design seed. Governs the site draws, the participant draws and the run order.
set.seed(20260827)

# The results directory is VERSIONED for the same reason as permutation_v3/ and
# full_v2/: a cell is skipped when its file exists, so a changed design or a
# changed set of recorded columns must not be written into an existing
# directory. Bump the suffix when the design, the arms, or the recorded meta
# columns change.
out_dir <- here::here("analysis/data/derived_data/timing_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The full run this experiment is validated against.
full_dir <- here::here("analysis/data/derived_data/full_v2")
full_timing_file <- file.path(full_dir, "timing.RDS")

# THE CORE COUNT IS PART OF THE MEASUREMENT. full_v2/timing.RDS records the
# core count the full run actually used, and a timing experiment run on a
# different number of cores does not predict it: the slope is per-target cost
# at that width, and the intercept is largely cluster start-up, which grows
# with the number of workers. (Note also that the case-study text currently
# says 12 cores while the recorded meta says 16 -- the file is the record.)
#
# When the two disagree the experiment is still worth running -- linearity in k
# is a property of seqwrap, not of the machine -- but the extrapolation below
# is then a comparison across machines, and it is reported as such rather than
# silently.
cores_full <- if (file.exists(full_timing_file)) {
  unique(readRDS(full_timing_file)$cores)
} else {
  NA_integer_
}

if (!is.na(cores_full) && !identical(as.integer(cores_full), as.integer(CORES))) {
  message(
    "NOTE: the full run used ", cores_full, " cores, this machine offers ",
    CORES, ". The linearity result is unaffected; the extrapolation to the ",
    "full run is a cross-machine comparison and is reported with a ",
    "perfect-scaling adjustment alongside the raw prediction."
  )
}


## Design parameters #######################################################

# Target counts for the size block. The lowest levels are there to pin the
# intercept (the per-call overhead), not because anyone would run seqwrap on
# 100 targets; the highest is what the slope is mostly estimated from.
ks <- if (QUICK) c(50, 100, 250) else c(100, 250, 500, 1000, 2000, 4000)

# Replicates per cell. Each replicate is an independent draw of sites (and of
# participants, in the samples block), so the spread between replicates is the
# site-composition variance plus machine noise -- which is the error term the
# linear fit at the end is read against.
reps <- if (QUICK) 2 else 3

# Fixed k for the two blocks that do not vary k.
k_complexity <- if (QUICK) 250 else 2000
k_samples <- if (QUICK) 250 else 1000

# Participants retained in the samples block. 8 is the full data set.
n_participants <- c(4, 6, 8)

# The anchor cell, re-run every `anchor_every` cells. Chosen as the cheapest
# arm at a middling k: it has to be run often, so it has to be cheap.
ANCHOR <- list(arm = "m", k = if (QUICK) 100 else 1000)
anchor_every <- 10


## Load data ###############################################################
#
# Identical to the head of methylation-case-full.R: the normalized gset is
# built once, upstream, by methylation-case-study-dataprep.R.

gset_file <- here::here(
  "analysis/data/derived_data/seaborne-gset-normalized.RDS"
)

if (!file.exists(gset_file)) {
  stop(
    "seaborne-gset-normalized.RDS is missing. Run ",
    "analysis/R/methylation-case-study-dataprep.R first.",
    call. = FALSE
  )
}

gset <- readRDS(gset_file)
metadata <- readRDS(
  here::here("analysis/data/derived_data/seaborne-metadata.RDS")
)

beta_vals <- getBeta(gset, offset = 100)
m_vals <- log2(beta_vals / (1 - beta_vals))

rm(gset)
gc()

betadf <- data.frame(beta_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

mdf <- data.frame(m_vals) |>
  tibble::rownames_to_column(var = "id") |>
  dplyr::select(id, starts_with("GSM"))

rm(beta_vals, m_vals)
gc()

# The two scales must be two representations of the same sites in the same
# order: the size block draws ONE set of row indices per replicate and uses it
# on both, so that the beta and M arms of a cell are timed on the same sites.
stopifnot(
  identical(betadf$id, mdf$id),
  identical(names(betadf), names(mdf))
)

n_sites <- nrow(betadf)
sample_cols <- setdiff(names(betadf), "id")

stopifnot(setequal(sample_cols, metadata$geo_accession))


## Model functions #########################################################
#
# eval_fun and the two summary functions below are kept IDENTICAL to
# methylation-case-full.R. That script is the source of truth; it cannot be
# sourced here because sourcing it fits the full data, so they are duplicated.
# If the full-run versions change, these have to change with them -- otherwise
# the timings are timings of different work and the extrapolation at the end of
# this script measures the divergence rather than the scaling.

eval_fun <- function(m) {
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

# The complexity ladder includes formulas with no participant random effect and
# with a second random effect, so the arm eval_fun above -- which reads
# VarCorr()$cond$participant -- does not apply to it. This is the version used
# by the complexity block, and it is used by EVERY rung of that ladder,
# including the rung that is the full-run model, so the ladder is internally
# comparable. It is deliberately not used elsewhere: the size and samples
# blocks are timed with the full run's own eval_fun.
eval_fun_timing <- function(m) {
  if (!inherits(m, "glmmTMB")) {
    stop("unsupported model class: ", paste(class(m), collapse = "/"))
  }

  vc <- glmmTMB::VarCorr(m)$cond

  data.frame(
    convergence = m$fit$convergence,
    singular = !m$sdr$pdHess,
    re_sd = if (length(vc)) sqrt(as.numeric(vc[[1]][1, 1])) else NA_real_,
    n = nrow(m$frame),
    reml = isTRUE(m$modelInfo$REML)
  )
}

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

  i <- if (is.null(sat)) NA_integer_ else match(out$term, rownames(sat))

  out$ddf <- if (is.null(sat)) NA_real_ else sat[i, "ddf"]
  out$p_satt <- if (is.null(sat)) NA_real_ else sat[i, "Pr(>|t|)"]

  out
}

summary_fun_wald_nosatt <- function(m) {
  cf <- coef(summary(m))$cond

  data.frame(
    term = rownames(cf),
    estimate = cf[, "Estimate"],
    std.error = cf[, "Std. Error"],
    p_wald = cf[, "Pr(>|z|)"],
    ddf = NA_real_,
    p_satt = NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


## Arms ####################################################################
#
# The four arms of the full run, verbatim. `data` is referenced by scale rather
# than stored, so the arm table does not hold two copies of a 230 MB frame.
arms <- list(
  beta = list(
    name = "beta",
    scale = "beta",
    summary_fun = summary_fun_wald_nosatt,
    arguments = list(
      formula = y ~ time + (1 | participant),
      family = glmmTMB::beta_family()
    )
  ),
  beta_reml = list(
    name = "beta_reml",
    scale = "beta",
    summary_fun = summary_fun_wald,
    arguments = list(
      formula = y ~ time + (1 | participant),
      family = glmmTMB::beta_family(),
      REML = TRUE
    )
  ),
  m = list(
    name = "m",
    scale = "M",
    summary_fun = summary_fun_wald_nosatt,
    arguments = list(
      formula = y ~ time + (1 | participant)
    )
  ),
  m_reml = list(
    name = "m_reml",
    scale = "M",
    summary_fun = summary_fun_wald,
    arguments = list(
      formula = y ~ time + (1 | participant),
      REML = TRUE
    )
  )
)

# The complexity ladder. All rungs are ML and are read with the plain Wald
# summary, so the ONLY thing that varies along the ladder is the model: the
# estimator and the reference distributions are held fixed. `age` is a
# between-participant covariate and `slide` the array batch; both are in the
# metadata. The ladder is run on both scales, since a beta likelihood and a
# Gaussian one need not respond to added structure in the same way.
ladder <- list(
  list(label = "fixed only", formula = y ~ time),
  list(label = "intercept + RE", formula = y ~ 1 + (1 | participant)),
  list(label = "time + RE", formula = y ~ time + (1 | participant)),
  list(label = "time + age + RE", formula = y ~ time + age + (1 | participant)),
  list(
    label = "time + 2 RE",
    formula = y ~ time + (1 | participant) + (1 | slide)
  )
)


## Sub-sample draws ########################################################
#
# Drawn once, here, and saved, so the design is a recorded object rather than
# something that has to be reproduced by re-running the seeds in order.

# Sites: one draw of max(k) row indices per replicate; the k levels are nested
# prefixes of that draw (see header).
k_max <- max(c(ks, k_complexity, k_samples, ANCHOR$k))
stopifnot(k_max <= n_sites)

site_draws <- lapply(seq_len(reps), function(r) {
  withr::with_seed(1000 + r, sample.int(n_sites, k_max))
})

# Participants: nested in the same way, so the 4-participant draw of a
# replicate is a subset of its 6-participant draw. A participant is dropped
# whole, with all their time points.
participants <- sort(unique(metadata$participant))

part_draws <- lapply(seq_len(reps), function(r) {
  ord <- withr::with_seed(2000 + r, sample(participants))
  lapply(
    stats::setNames(n_participants, n_participants),
    function(np) sort(ord[seq_len(np)])
  )
})

saveRDS(
  list(
    site_draws = site_draws,
    part_draws = part_draws,
    site_ids = betadf$id[unique(unlist(site_draws))]
  ),
  file.path(out_dir, "draws.RDS")
)


## The design table ########################################################
#
# One row per seqwrap() call. `cell` is the file key; it is built from the
# design fields rather than from the row number, so that adding a block later
# does not renumber -- and therefore invalidate -- the cells already on disk.

design_size <- expand.grid(
  arm = names(arms),
  k = ks,
  rep = seq_len(reps),
  stringsAsFactors = FALSE
) |>
  mutate(block = "size", ladder = NA_integer_, n_part = max(n_participants))

design_complexity <- expand.grid(
  scale = c("beta", "M"),
  ladder = seq_along(ladder),
  rep = seq_len(reps),
  stringsAsFactors = FALSE
) |>
  mutate(
    block = "complexity",
    arm = if_else(scale == "beta", "beta", "m"),
    k = k_complexity,
    n_part = max(n_participants)
  ) |>
  select(-scale)

design_samples <- expand.grid(
  arm = names(arms),
  n_part = n_participants,
  rep = seq_len(reps),
  stringsAsFactors = FALSE
) |>
  mutate(block = "samples", ladder = NA_integer_, k = k_samples)

design <- bind_rows(design_size, design_complexity, design_samples) |>
  mutate(
    cell = sprintf(
      "%s-%s-k%05d-n%d-l%s-r%d",
      block,
      arm,
      k,
      n_part,
      if_else(is.na(ladder), "NA", sprintf("%d", ladder)),
      rep
    )
  )

stopifnot(!anyDuplicated(design$cell))

# Run order: randomised over the whole design, then the anchor spliced in at
# regular positions. Randomisation is what keeps a machine that slows down over
# two hours from looking like a design effect; the anchor is what detects that
# it did.
design <- design[sample.int(nrow(design)), ]

# Every anchor is the SAME cell: same arm, same k, same sites (replicate 1),
# same participants. That is the point of it -- anything that differs between
# two anchors is the machine, not the design, so the replicate index is fixed
# here rather than following the anchor number.
anchor_row <- function(i) {
  data.frame(
    block = "anchor",
    arm = ANCHOR$arm,
    k = ANCHOR$k,
    n_part = max(n_participants),
    ladder = NA_integer_,
    rep = 1L,
    cell = sprintf("anchor-%03d", i),
    stringsAsFactors = FALSE
  )
}

run_order <- list()
n_anchor <- 0L
for (i in seq_len(nrow(design))) {
  if ((i - 1L) %% anchor_every == 0L) {
    n_anchor <- n_anchor + 1L
    run_order[[length(run_order) + 1L]] <- anchor_row(n_anchor)
  }
  run_order[[length(run_order) + 1L]] <- design[i, , drop = FALSE]
}
run_order <- bind_rows(run_order)

saveRDS(run_order, file.path(out_dir, "design.RDS"))

message(
  "design: ",
  nrow(run_order),
  " cells (",
  n_anchor,
  " anchors) on ",
  CORES,
  " cores"
)


## Running one cell ########################################################

# Everything a cell needs that is not in the design row: the data frame for the
# scale, the metadata rows, the model arguments and the summary/eval functions.
cell_spec <- function(row) {
  arm <- arms[[row$arm]]

  keep_part <- part_draws[[row$rep]][[as.character(row$n_part)]]
  keep_meta <- metadata[metadata$participant %in% keep_part, , drop = FALSE]

  idx <- site_draws[[row$rep]][seq_len(row$k)]

  dat <- if (arm$scale == "beta") betadf else mdf
  # Columns as well as rows: seqwrap matches the response columns to the
  # metadata by `samplename`, so an array that is not in the metadata of this
  # cell must not be in its data either.
  dat <- dat[idx, c("id", keep_meta$geo_accession), drop = FALSE]

  if (row$block == "complexity") {
    rung <- ladder[[row$ladder]]
    arguments <- arm$arguments
    arguments$formula <- rung$formula
    arguments$REML <- NULL
    list(
      data = dat,
      metadata = keep_meta,
      arguments = arguments,
      summary_fun = summary_fun_wald_nosatt,
      eval_fun = eval_fun_timing,
      model_label = rung$label
    )
  } else {
    list(
      data = dat,
      metadata = keep_meta,
      arguments = arm$arguments,
      summary_fun = arm$summary_fun,
      eval_fun = eval_fun,
      model_label = arm$name
    )
  }
}

# The number of parameters the optimiser actually works on, measured rather
# than counted off the formula: under REML the fixed effects are integrated
# out, so `np_opt` is smaller than the number of coefficients reported, and
# that is part of why the REML arms are not slower for the reason one might
# assume. Fitted OUTSIDE the timed region, on one target, so it costs the
# experiment one fit per cell and contaminates nothing.
probe_np <- function(spec) {
  fit <- try(
    {
      d <- data.frame(
        y = as.numeric(spec$data[1, spec$metadata$geo_accession]),
        spec$metadata
      )
      do.call(glmmTMB::glmmTMB, c(list(data = d), spec$arguments))
    },
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(data.frame(np_opt = NA_integer_, np_fixed = NA_integer_))
  }

  data.frame(
    np_opt = length(fit$fit$par),
    np_fixed = length(glmmTMB::fixef(fit)$cond)
  )
}

run_cell <- function(row) {
  spec <- cell_spec(row)

  # Composition is outside seqwrap()'s own clock, so it is timed separately
  # rather than left unaccounted.
  t_compose <- system.time({
    container <- seqwrap_compose(
      modelfun = glmmTMB::glmmTMB,
      data = spec$data,
      metadata = spec$metadata,
      samplename = "geo_accession",
      eval_fun = spec$eval_fun,
      summary_fun = spec$summary_fun,
      arguments = spec$arguments
    )
  })

  np <- probe_np(spec)

  # The timed region. Arguments other than the data are those of the full run.
  results <- seqwrap(
    container,
    return_models = FALSE,
    cores = CORES,
    verbose = FALSE
  )

  summarised <- seqwrap_summarise(results, verbose = FALSE)

  ok <- !is.null(summarised$evaluations) && NROW(summarised$evaluations)

  meta <- data.frame(
    cell = row$cell,
    block = row$block,
    arm = row$arm,
    model_label = spec$model_label,
    ladder = row$ladder,
    scale = if (row$arm %in% c("beta", "beta_reml")) "beta" else "M",
    reml = isTRUE(spec$arguments$REML),
    satterthwaite = identical(spec$summary_fun, summary_fun_wald),
    rep = row$rep,
    cores = CORES,
    # Recorded because it is the parallel granularity and it is a function of
    # k: if the timings were ever non-linear in k this is the first place to
    # look.
    chunk_size = max(1L, min(2000L, ceiling(row$k / (CORES * 4L)))),
    k = results@k,
    n = results@n,
    n_part = row$n_part,
    np_opt = np$np_opt,
    np_fixed = np$np_fixed,
    compose = t_compose[["elapsed"]],
    user = results@elapsed_time[["user.self"]],
    system = results@elapsed_time[["sys.self"]],
    elapsed = results@elapsed_time[["elapsed"]],
    # A cell in which many fits failed is a cell that was cheap for the wrong
    # reason, so the failure rate travels with the timing.
    n_evaluated = if (ok) NROW(summarised$evaluations) else 0L,
    n_converged = if (ok) {
      sum(summarised$evaluations$convergence == 0, na.rm = TRUE)
    } else {
      0L
    },
    n_singular = if (ok) {
      sum(summarised$evaluations$singular, na.rm = TRUE)
    } else {
      0L
    },
    n_conditions = NROW(summarised$errors),
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )

  # The design row asked for k targets; seqwrap reports what it fitted.
  stopifnot(meta$k == row$k, meta$n == nrow(spec$metadata))

  rm(container, results, summarised)
  gc()

  meta
}


## Warm-up #################################################################
#
# The first seqwrap() call of a session pays for things no later call pays for
# again: package load on the workers, byte-compilation, a cold file cache. It
# is run and discarded so that whichever cell the randomisation happens to put
# first is not the one that carries that cost.
if (!file.exists(file.path(out_dir, "warmup.RDS"))) {
  message("warm-up run")
  warm <- run_cell(
    data.frame(
      block = "warmup",
      arm = "m",
      k = min(ks),
      n_part = max(n_participants),
      ladder = NA_integer_,
      rep = 1L,
      cell = "warmup",
      stringsAsFactors = FALSE
    )
  )
  saveRDS(warm, file.path(out_dir, "warmup.RDS"))
}


## Run the experiment ######################################################
#
# One file per cell, skipped when it exists.
for (i in seq_len(nrow(run_order))) {
  row <- run_order[i, , drop = FALSE]
  cell_file <- file.path(out_dir, sprintf("cell-%s.RDS", row$cell))

  if (file.exists(cell_file)) next

  message(
    sprintf(
      "[%3d/%3d] %s  k=%d  n=%d",
      i,
      nrow(run_order),
      row$cell,
      row$k,
      row$n_part
    )
  )

  meta <- run_cell(row)
  saveRDS(meta, cell_file)

  message(sprintf("           %.1f s elapsed", meta$elapsed))
}


## Collect #################################################################

cell_files <- list.files(
  out_dir,
  pattern = "^cell-.*\\.RDS$",
  full.names = TRUE
)

timing <- bind_rows(lapply(cell_files, readRDS)) |>
  mutate(
    estimator = if_else(reml, "REML", "ML"),
    # The quantity the whole experiment is about.
    per_target = elapsed / k,
    # Core-seconds per target: what a single core would have spent. Not the
    # same claim as `per_target`, and the gap between them is the parallel
    # efficiency.
    core_per_target = elapsed * cores / k
  )

saveRDS(timing, file.path(out_dir, "timing-cells.RDS"))


## Analysis ################################################################

# 1. Is the cost linear in k?
#
# Two readings of the same size block. The linear fit gives the two numbers
# that are the practical answer -- a fixed cost per call (the intercept, in
# seconds) and a marginal cost per target (the slope) -- but it ASSUMES the
# linearity it is being used to describe. The log-log fit is the test:
# elapsed ~ k^b, and linearity is b = 1. A super-linear cost (b > 1) is what
# would make an experiment on 4 000 targets useless for predicting a run on
# 770 441.
size <- filter(timing, block == "size")

# A design with only two k levels and one replicate leaves no residual degrees
# of freedom, and confint() then returns NaN with a warning rather than
# failing. That is the QUICK design, not the experiment, but the guard is here
# so that a smoke test reports "no interval" instead of a NaN that reads like a
# result.
ci_or_na <- function(m, which) {
  if (stats::df.residual(m) < 1) {
    return(NA_real_)
  }
  confint(m)["log(k)", which]
}

scaling <- size |>
  group_by(arm, scale, estimator) |>
  summarise(
    lin = list(lm(elapsed ~ k)),
    loglog = list(lm(log(elapsed) ~ log(k))),
    .groups = "drop"
  ) |>
  mutate(
    overhead_s = vapply(lin, function(m) coef(m)[[1]], numeric(1)),
    per_target_s = vapply(lin, function(m) coef(m)[[2]], numeric(1)),
    r2 = vapply(lin, function(m) summary(m)$r.squared, numeric(1)),
    exponent = vapply(loglog, function(m) coef(m)[[2]], numeric(1)),
    exponent_lo = vapply(loglog, ci_or_na, numeric(1), which = 1),
    exponent_hi = vapply(loglog, ci_or_na, numeric(1), which = 2),
    # The reading that matters: does the confidence interval on the exponent
    # contain 1?
    linear = exponent_lo <= 1 & exponent_hi >= 1
  )


# 2. Does the small-scale fit predict the full run?
#
# This is the check that decides whether the experiment says anything about the
# 10-hour run. The prediction is made from cells of at most 4 000 targets and
# is compared with the observed full-run timing at 770 441 -- an extrapolation
# of more than two orders of magnitude, which nothing in the fit is protecting.
#
# `predicted_h` is the honest prediction: what a run of k_full targets would
# cost ON THIS MACHINE, at this core count. `predicted_h_adj` rescales it to
# the core count the full run used, and it is only meaningful when the two
# differ; the rescaling assumes perfect scaling in cores, which is false, so it
# is a bound on the comparison rather than a prediction. When the core counts
# agree the two columns are identical and only the first is worth reading.
if (file.exists(full_timing_file)) {
  full_timing <- readRDS(full_timing_file)

  extrapolation <- scaling |>
    select(arm, overhead_s, per_target_s) |>
    left_join(
      transmute(
        full_timing,
        arm = model,
        observed_h = hours,
        k_full = k,
        cores_full = cores
      ),
      by = "arm"
    ) |>
    mutate(
      cores_here = CORES,
      predicted_h = (overhead_s + per_target_s * k_full) / 3600,
      predicted_h_adj = predicted_h * cores_here / cores_full,
      ratio = predicted_h_adj / observed_h,
      error_pct = 100 * (predicted_h_adj - observed_h) / observed_h
    )
} else {
  message("full_v2/timing.RDS not found -- skipping the extrapolation check")
  extrapolation <- NULL
}


# 3. Does the cost track model complexity?
#
# Per-target cost along the ladder, at fixed k and fixed sample size, against
# the number of parameters the optimiser worked on. `np_opt` is a crude summary
# of a model -- it does not distinguish a random effect from a fixed one, and
# the random-effect terms are what actually drive the cost of a Laplace
# approximation -- so the ladder labels are kept alongside it rather than being
# reduced to it.
complexity <- timing |>
  filter(block == "complexity") |>
  group_by(scale, ladder, model_label, np_opt, np_fixed) |>
  summarise(
    per_target_ms = 1000 * mean(per_target),
    sd_ms = 1000 * sd(per_target),
    converged = sum(n_converged) / sum(n_evaluated),
    .groups = "drop"
  ) |>
  arrange(scale, ladder)


# 4. Does the cost track the number of arrays?
samples <- timing |>
  filter(block == "samples") |>
  group_by(arm, scale, estimator, n_part, n) |>
  summarise(
    per_target_ms = 1000 * mean(per_target),
    sd_ms = 1000 * sd(per_target),
    .groups = "drop"
  ) |>
  arrange(arm, n_part)


# 5. Did the machine hold still?
#
# The anchor cells in run order. A slope here is machine drift, and it is the
# thing that would otherwise be read as an effect of whatever the design
# happened to schedule late.
anchor <- timing |>
  filter(block == "anchor") |>
  arrange(timestamp) |>
  mutate(order = row_number())

anchor_drift <- if (nrow(anchor) > 2) {
  fit <- lm(elapsed ~ order, data = anchor)
  data.frame(
    n_anchor = nrow(anchor),
    mean_s = mean(anchor$elapsed),
    cv_pct = 100 * sd(anchor$elapsed) / mean(anchor$elapsed),
    slope_s_per_cell = coef(fit)[[2]],
    p_value = summary(fit)$coefficients["order", 4]
  )
} else {
  NULL
}


saveRDS(
  list(
    scaling = mutate(scaling, lin = NULL, loglog = NULL),
    extrapolation = extrapolation,
    complexity = complexity,
    samples = samples,
    anchor = anchor_drift,
    session = list(
      r_version = R.version.string,
      platform = R.version$platform,
      cores = CORES,
      seqwrap = as.character(utils::packageVersion("seqwrap")),
      glmmTMB = as.character(utils::packageVersion("glmmTMB"))
    )
  ),
  file.path(out_dir, "timing-scaling.RDS")
)


## Figure ##################################################################
#
# Left: elapsed against k on linear axes with the fitted line -- the reading a
# user cares about. Right: the same on log axes with slope-1 reference lines --
# the reading that says the left panel is entitled to be a straight line.

p_linear <- ggplot(size, aes(k, elapsed, colour = arm)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  geom_point(alpha = 0.7) +
  labs(x = "Targets (k)", y = "Elapsed (s)", colour = "Arm") +
  theme_minimal()

p_log <- ggplot(size, aes(k, elapsed, colour = arm)) +
  geom_abline(
    slope = 1,
    intercept = seq(-8, 0, by = 1),
    colour = "grey85",
    linewidth = 0.3
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Targets (k, log)",
    y = "Elapsed (s, log)",
    colour = "Arm",
    caption = "grey lines: slope 1 (cost proportional to k)"
  ) +
  theme_minimal()

fig <- cowplot::plot_grid(
  p_linear + theme(legend.position = "none"),
  p_log,
  nrow = 1,
  rel_widths = c(1, 1.25)
)

ggsave(
  file.path(out_dir, "timing-scaling.png"),
  fig,
  width = 10,
  height = 4,
  dpi = 300
)


## Report ##################################################################

print(as.data.frame(select(
  scaling,
  arm,
  overhead_s,
  per_target_s,
  r2,
  exponent,
  exponent_lo,
  exponent_hi,
  linear
)))

if (!is.null(extrapolation)) {
  cat("\nExtrapolation to the full run:\n")
  print(as.data.frame(select(
    extrapolation,
    arm,
    k_full,
    cores_here,
    cores_full,
    predicted_h,
    predicted_h_adj,
    observed_h,
    error_pct
  )))
}

cat("\nComplexity ladder (per-target ms):\n")
print(as.data.frame(complexity))

cat("\nSample size (per-target ms):\n")
print(as.data.frame(samples))

if (!is.null(anchor_drift)) {
  cat("\nAnchor cells (machine drift):\n")
  print(anchor_drift)
}
