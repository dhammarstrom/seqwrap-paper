# Post-fitting shrinkage of the methylation effect estimates ##################
#
# An example of what the full run supports once it exists. seqwrap returns an
# estimate and a standard error per target, which is exactly the input that
# post-fitting machinery takes, so the case study can end by declaring a set of
# differentially methylated positions rather than by tabulating p-values.
#
# The method is adaptive shrinkage (Stephens 2017, Biostatistics; the ashr
# package). It assumes
#
#     (beta_hat_j - beta_j) / se_j  ~  t(df_j)
#
# and is normally applied on the hope that this holds. Here it does not have to
# be hoped for: the permutation study (analysis/R/methylation-error-permutation.R)
# established that the REML standard errors read against a t with Satterthwaite
# denominator df are the calibrated ones, and those are the two columns fed in
# below. That is the reason this step is applied to the REML arms and to no
# others, and it is why nothing further is validated here.
#
# The procedure itself is one call per (arm, contrast):
#
#   ash(estimate, std.error, mixcompdist = "halfuniform", lik = lik_t(df = ddf))
#
# with the number of declarations read off as sum(get_lfsr(fit) <= 0.05).
# Four points about that call, each marked again where it is made:
#
# (1) DEGREES OF FREEDOM. ash() takes a scalar `df` and stops on a vector
#     ("Only one value can be specified for df."). The Satterthwaite ddf is per
#     site, so it goes in through lik_t(), which does take a vector -- it
#     carries `const = length(unique(df)) == 1` precisely to note when the df
#     varies. This is the only non-default argument that matters.
#
# (2) SCALE. Both arms are put in M units first. The beta arm is a logit-link
#     model and logit(beta) = M * log(2), so its estimate AND its standard
#     error are divided by log(2) -- the same conversion `est_M` uses in the
#     case study. Every t statistic is unchanged by it. Note that a posterior
#     mean cannot then be pushed back through plogis() and called a posterior
#     mean of a difference in beta: the transform is nonlinear.
#
# (3) ONE FIT PER CONTRAST. The prior is the distribution of true effects
#     across sites, and the acute response and the reloading response have no
#     reason to share one.
#
# (4) DECLARATIONS COME FROM THE LFSR, not from the s-value. The lfsr is the
#     posterior probability of getting a site's sign wrong and is the per-site
#     quantity. The s-value is its cumulative average over the ranked list, so
#     it saturates near one on a contrast with a large non-null fraction and
#     is not a per-site threshold. Both are stored; `n_lfsr_05` is the count
#     the case study reports.
#
# ONE THING TO READ BEFORE THE COUNTS, because on two of these four contrasts
# `n_lfsr_05` does not measure what it appears to.
#
# ash assumes the true effects are unimodal about zero. 2_acute and 3_loading
# are not centred at zero -- the mean estimate is -0.023 and -0.034 M units,
# against -0.003 at 4_unloading and +0.011 at 5_reloading -- and a prior with
# its mode pinned at zero can only explain a bulk sitting off zero by putting
# all of its non-null mass on one side. That is exactly what it does: on both
# contrasts the fitted prior is 100% negative and every posterior mean comes
# back negative.
#
# A one-sided prior puts a FLOOR under the lfsr. For a site with no evidence of
# its own the posterior is the prior, so P(beta >= 0) = pi0 and P(beta <= 0) = 1,
# giving lfsr = pi0 exactly. Measured: pi0 = 0.077 and lfsr = 0.079 at 2_acute,
# pi0 = 0.042 and lfsr = 0.043 at 3_loading. Thresholding at 0.05 is then a
# switch rather than a test -- below pi0 it declares almost nothing, above pi0
# almost everything -- and those two contrasts happen to fall on opposite sides
# of it, which is the whole of the difference between declaring 4% of sites and
# declaring 84%. Neither number is a count of sites with per-site evidence.
#
# `prior_oneside` and `lfsr_uninformative` are in the summary table for this
# reason. When `prior_oneside` is at or near 1 and `lfsr_uninformative` sits
# near the threshold, the count is reporting 1 - pi0 and nothing else. On
# 4_unloading and 5_reloading the prior is genuinely two-sided (a 0.71 and a
# 0.63 majority side), the uninformative lfsr is 0.80 and 0.69, and the counts
# mean what they say.
#
# The offset that causes this is not something this script can resolve. A shift
# shared by every probe on the array is either a real global change in
# methylation at that timepoint or residual technical variation between
# timepoints, and the two enter the effect estimates identically. It belongs
# upstream, with the normalization. Re-centring the prior instead
# (mode = "estimate") is NOT a fix: it drives pi0 to ~1 and makes the lfsr
# meaningless in the other direction, because once the prior mode sits at -0.03
# every posterior is confidently negative.
#
# Outputs, under the versioned directory set by `out_dir`:
#
#   shrink-<arm>-<term>.RDS   per-site posterior summaries, plus the fitted
#                             prior, for one arm and contrast.
#   shrinkage-summary.RDS     one row per arm x contrast: pi0, the counts, and
#                             how far the shrinkage moved the estimates.
#
# A file that already exists is skipped, so an interrupted run resumes at the
# contrast it stopped in. Bump the directory suffix when the arms, the prior or
# the output schema change, so a new format cannot be mixed into an old
# directory.
#
# Cost: about three and a half minutes per contrast at this number of sites,
# so under half an hour for the eight fits, plus the arm reads.


## 00. Packages #############################################################

library(dplyr)

# ashr is not yet declared by the compendium: it appears in renv.lock only as a
# requirement of DESeq2. `etrunct` is what lik_t() uses for the truncated-t
# moments and comes in with it.
#
#   renv::install("ashr")   then   renv::snapshot()
#
# and add ashr to Imports in DESCRIPTION.
for (pkg in c("ashr", "etrunct")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      pkg,
      " is not installed. Run renv::install(\"ashr\"), then renv::snapshot().",
      call. = FALSE
    )
  }
}


## 01. Configuration ########################################################

full_dir <- here::here("analysis/data/derived_data/full_v2")
out_dir <- here::here("analysis/data/derived_data/shrinkage_v1")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(full_dir))

# The REML arms only, for the reason given in the header. `divisor` puts the
# arm in M units, point (2).
shrink_arms <- list(
  list(name = "m_reml", scale = "M", divisor = 1),
  list(name = "beta_reml", scale = "beta", divisor = log(2))
)

# 1_baseline is the reference level; the intercept is not a contrast.
full_terms <- c(
  "time2_acute",
  "time3_loading",
  "time4_unloading",
  "time5_reloading"
)

# "halfuniform" keeps the unimodal assumption but lets the two sides of the
# prior carry their own weight, which a training intervention has no reason to
# balance. "normal" and "halfnormal" are unavailable with a t likelihood.
ash_mixcompdist <- "halfuniform"

# The lfsr threshold the counts are read at.
lfsr_threshold <- 0.05

# A site enters only if its fit converged, was not singular, and returned a
# usable estimate, standard error and ddf. The ddf floor is 2 because that is
# where the t stops having a variance: below it ash still returns a posterior
# mean but the posterior SD is NaN. It costs nothing here -- the smallest ddf
# across the full run is 13.1 -- but a boundary fit should not be reported with
# an uncertainty that is not defined.
ddf_min <- 2


## 02. Shrink each arm and contrast #########################################
#
# The arms are read one at a time and dropped again: an arm file holds 3.85
# million summary rows and two of them alive at once is what pushes this into
# swap on a laptop.

for (arm in shrink_arms) {
  arm_file <- file.path(full_dir, sprintf("%s.RDS", arm$name))

  out_files <- file.path(
    out_dir,
    sprintf("shrink-%s-%s.RDS", arm$name, full_terms)
  )

  if (all(file.exists(out_files))) {
    message("skipping arm ", arm$name, " -- all contrasts already shrunk")
    next
  }

  if (!file.exists(arm_file)) {
    stop(
      "missing ",
      basename(arm_file),
      ". Run analysis/R/methylation-case-full.R first.",
      call. = FALSE
    )
  }

  message("reading arm ", arm$name)
  full <- readRDS(arm_file)

  stopifnot(
    all(
      c("target", "term", "estimate", "std.error", "ddf") %in%
        names(full$summaries)
    ),
    all(c("target", "convergence", "singular") %in% names(full$evaluations)),
    setequal(setdiff(unique(full$summaries$term), "(Intercept)"), full_terms),
    # `reml` is read off the fit object rather than the arm name, so this is
    # the check that the arm is the one it is named after.
    isTRUE(unique(full$evaluations$reml))
  )

  for (trm in full_terms) {
    out_file <- file.path(out_dir, sprintf("shrink-%s-%s.RDS", arm$name, trm))

    if (file.exists(out_file)) {
      message("  skipping ", arm$name, " / ", trm)
      next
    }

    dat <- full$summaries |>
      dplyr::filter(term == trm) |>
      dplyr::inner_join(
        dplyr::select(full$evaluations, target, convergence, singular),
        by = "target"
      ) |>
      dplyr::filter(
        convergence == 0,
        !singular,
        is.finite(estimate),
        is.finite(std.error),
        std.error > 0,
        is.finite(ddf),
        ddf >= ddf_min
      ) |>
      dplyr::transmute(
        target = target,
        # M units, point (2). Estimate and standard error are divided by the
        # same constant, so the t statistics are untouched.
        estimate_M = estimate / arm$divisor,
        se_M = std.error / arm$divisor,
        ddf = ddf
      )

    stopifnot(nrow(dat) > 0, !anyDuplicated(dat$target))

    message("  shrinking ", arm$name, " / ", trm, " (", nrow(dat), " sites)")

    # The procedure. Point (1): the per-site ddf enters through lik_t(),
    # because ash(df = ) takes a scalar only.
    elapsed <- system.time(
      fit <- ashr::ash(
        dat$estimate_M,
        dat$se_M,
        mixcompdist = ash_mixcompdist,
        lik = ashr::lik_t(df = dat$ddf)
      )
    )

    out <- dat |>
      dplyr::mutate(
        post_mean = ashr::get_pm(fit),
        post_sd = ashr::get_psd(fit),
        lfsr = ashr::get_lfsr(fit),
        svalue = ashr::get_svalue(fit)
      )

    # A NaN reaching a saved file would surface much later, in a figure.
    stopifnot(
      nrow(out) == nrow(dat),
      all(is.finite(out$post_mean)),
      all(is.finite(out$lfsr))
    )

    saveRDS(
      list(
        model = arm$name,
        term = trm,
        scale = arm$scale,
        dat = out,
        pi0 = ashr::get_pi0(fit),
        fitted_g = ashr::get_fitted_g(fit),
        elapsed = elapsed
      ),
      out_file
    )

    rm(fit, out, dat)
    gc()
  }

  rm(full)
  gc()
}


## 03. The table the case study reports #####################################
#
# One row per arm x contrast. `n_lfsr_05` is the count of declared positions;
# `mean_est` is the genome-wide offset that has to be read alongside it, see
# the header.

summary_file <- file.path(out_dir, "shrinkage-summary.RDS")

shrink_files <- list.files(
  out_dir,
  pattern = "^shrink-.*[.]RDS$",
  full.names = TRUE
)

if (!file.exists(summary_file) && length(shrink_files) > 0) {
  shrinkage_summary <- dplyr::bind_rows(lapply(shrink_files, function(f) {
    s <- readRDS(f)
    d <- s$dat
    declared <- d$lfsr <= lfsr_threshold

    out <- data.frame(
      model = s$model,
      scale = s$scale,
      term = s$term,
      n = nrow(d),
      # The offset, carried next to the count it explains.
      mean_est = mean(d$estimate_M),
      pi0 = s$pi0,
      # The two columns that say whether n_lfsr_05 is measuring evidence or
      # only the prior, see the header.
      #
      # prior_oneside: the share of the non-null prior mass falling on
      # whichever side is heavier. At 1 the prior has concluded that every
      # real effect has the same sign.
      # lfsr_uninformative: what a site with no evidence of its own is given.
      # Under a fully one-sided prior this equals pi0 exactly, and it is the
      # floor of the whole lfsr distribution.
      prior_oneside = {
        g <- s$fitted_g
        mid <- (g$a + g$b) / 2
        neg <- sum(g$pi[mid < 0])
        pos <- sum(g$pi[mid > 0])
        if (neg + pos == 0) NA_real_ else max(neg, pos) / (neg + pos)
      },
      lfsr_uninformative = {
        flat <- abs(d$estimate_M / d$se_M) < 0.1
        if (any(flat)) stats::median(d$lfsr[flat]) else NA_real_
      },
      n_lfsr_05 = sum(declared),
      # What the shrinkage did to what is reported: the magnitude of the
      # declared effects before and after, in M units.
      median_abs_est_declared = if (any(declared)) {
        stats::median(abs(d$estimate_M[declared]))
      } else {
        NA_real_
      },
      median_abs_pm_declared = if (any(declared)) {
        stats::median(abs(d$post_mean[declared]))
      } else {
        NA_real_
      },
      elapsed = s$elapsed[["elapsed"]],
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    rm(s, d)
    gc()

    out
  }))

  shrinkage_summary <- dplyr::arrange(shrinkage_summary, model, term)

  saveRDS(shrinkage_summary, summary_file)

  print(shrinkage_summary[, c(
    "model",
    "term",
    "mean_est",
    "pi0",
    "prior_oneside",
    "lfsr_uninformative",
    "n_lfsr_05",
    "median_abs_pm_declared"
  )])
} else if (file.exists(summary_file)) {
  message("skipping the summary -- ", basename(summary_file), " exists")
}
