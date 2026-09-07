# Methylation experiment analysis #############################################
#
#
# 00. Loading package and dependencies
# 01. Beta categories used for the permutation test
# 02. Type I error rates per beta-category and p-valu distributions
# 02.B Upper panel
# 03. Standard errors per category and inferential procedure
# 04. Power (only in the REML-Satterthwaite procedure)
# 05. Results from full run
# 06. Timing experiment
# 07. Combine full figure
##############################################################################

# 00. Loading package and dependencies #######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(maxprobes) # For filtering cross-reactive probes
library(seqwrap)


# Load color scale
source(here::here("analysis/figures/figure-opts.R"))


# Label functions

# Labels shared by the table and the figures below, matching those used for
# the type I error rate.
lab_model <- function(x) {
  dplyr::recode(
    as.character(x),
    beta = "&beta;<sub>ML</sub>",
    beta_reml = "&beta;<sub>REML</sub>",
    m = "M<sub>ML</sub>",
    m_reml = "M<sub>REML</sub>"
  )
}

lab_type <- function(x) {
  factor(
    dplyr::recode(
      as.character(x),
      p_wald = "Wald p-values",
      p_satt = "Satterthwaite DDF",
      LRT = "Omnibus LRT"
    ),
    levels = c("Wald p-values", "Satterthwaite DDF", "Omnibus LRT")
  )
}

lab_betacat <- function(x) {
  factor(
    dplyr::case_when(
      x == "xlow" ~ "0 - 0.2",
      x == "low" ~ "0.2 - 0.4",
      x == "mid" ~ "0.4 - 0.6",
      x == "high" ~ "0.6 - 0.8",
      x == "xhigh" ~ "0.8 - 1"
    ),
    levels = c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1")
  )
}


axis_title_size <- 8
axis_text_size <- 8


# 01. Beta categories used for the permutation test  #########################

# list files in permutation results folder
# Versioned to match the script: permutation/ holds the earlier two-arm run
# (beta, m) with the old evaluation column names, permutation_v2/ the
# four-arm run (beta, beta_reml, m, m_reml_satt).
out_dir <- here::here("analysis/data/derived_data/permutation_v4")
perm_files <- list.files(
  out_dir,
  pattern = "^perm_\\d+\\.RDS$",
  full.names = TRUE
)

if (length(perm_files) == 0)
  source(
    here::here("analysis/R/methylation-error-permutation.R")
  )

n_perm <- length(perm_files)


### Reading the results back
spike_time <- "3_loading"

perm_files <- list.files(
  out_dir,
  pattern = "^perm_\\d+\\.RDS$",
  full.names = TRUE
)
perm <- lapply(perm_files, readRDS)

perm_summaries <- dplyr::bind_rows(lapply(perm, `[[`, "summaries"))
perm_designs <- dplyr::bind_rows(lapply(perm, `[[`, "design"))
perm_evals <- dplyr::bind_rows(lapply(perm, `[[`, "evaluations"))


# Which coefficients with the simulated effect.
#
# The added effect shifts only the samples labelled `spike_time` under the
# permuted labels. `spike_time` is set to 3_loading and
# 1_baseline is the reference level.

# alt_terms also accepts and omnibus test (not used here)
# Instead, the omnibus test is saved per target in the eval slot.
alt_terms <- c(paste0("time_permute", spike_time))

# Attach the per-iteration effects (delta, sign, stratum) to every test.
# is_alt is ground truth for a (site, term) pair, not for a site.
perm_dat <- perm_summaries |>
  dplyr::inner_join(perm_designs, by = c("iter", "target" = "id")) |>
  # Carry the per-fit convergence code across from the evaluation slot so
  # that failed fits can be dropped from the Wald error rates too. This
  # matters because the arms can fail at different rates, which would
  # quietly change which sites each arm is scored on.
  dplyr::inner_join(
    perm_evals |>
      dplyr::select(iter, target, model, convergence),
    by = c("iter", "target", "model")
  ) |>
  # delta is signed, so the test must be != 0 -- `delta > 0` would classify
  # every downward-spiked site as a null and contaminate the reference.
  dplyr::mutate(is_alt = delta != 0 & term %in% alt_terms)

stopifnot(any(perm_dat$is_alt))


# This figure extracts the &beta; values and shows the strata
# used in  the permutation study.

# Extract the beta values  ##############

# Reloading gset to make code chunk independent of above
gset <- readRDS(
  here::here("analysis/data/derived_data/seaborne-gset-quantile.RDS")
)

# Beta derived from M, not getBeta(offset = 100): the GenomicRatioSet holds no
# intensities and would ignore the offset silently. M is materialised once --
# it is a 10^2 MB matrix.
Mmat <- minfi::getM(gset)
B <- 2^Mmat / (1 + 2^Mmat)
rm(Mmat)

# Pre-compute densities and save in data frame
d <- density(rowMeans(B), from = 0, to = 1)
dd <- data.frame(x = d$x, y = d$y)


# Clwean up
rm(gset)
rm(B)
gc()


fig1 <- dd |>
  ggplot(aes(x, y)) +

  # Fill areas under
  geom_area(data = subset(dd, x <= 0.2), fill = colors[1], alpha = .6) +
  geom_area(
    data = subset(dd, x >= 0.2 & x <= 0.4),
    fill = colors[2],
    alpha = .6
  ) +
  geom_area(
    data = subset(dd, x >= 0.4 & x <= 0.6),
    fill = colors[3],
    alpha = .6
  ) +
  geom_area(
    data = subset(dd, x >= 0.6 & x <= 0.8),
    fill = colors[4],
    alpha = .6
  ) +
  geom_area(
    data = subset(dd, x >= 0.8 & x <= 1),
    fill = colors[5],
    alpha = .6
  ) +

  geom_line() +

  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    plot.subtitle = ggtext::element_markdown(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = colors)
  ) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    breaks = seq(from = 0, to = 1, by = 0.2),
    expand = c(0, 0)
  ) +
  labs(subtitle = "&beta;-values")


# 02. Type I error rates per beta-category and p-valu distributions ###########

# First we extract null targets among the wald-style parameters

null_dat <- perm_dat |>
  # Keep null effects in non affected sites
  filter(!(delta != 0 & term == "time_permute3_loading")) |>
  # filter NA in p-values
  filter(!is.na(p_satt))

# Alternative, keep only all nulls
# filter(delta == 0)

# The LR-test was removed in the latest permutation runs
# to simplify presentation. LR-test are not recommended for small sample
# glmm (see Bolker 2009).

# Collect LRT test statistics
# lrt_null <- perm_evals |>
#   filter(convergence == 0) |>
#   dplyr::inner_join(perm_designs, by = c("iter", "target" = "id")) |>
#   filter(delta == 0) |>
#
#   # fail safe if mssing p-vals from larger
#   # simulation
#   filter(!is.na(pval)) |>
#   mutate(
#     Estimator = "ML",
#     Model = if_else(model == "beta", "&beta;<sub>ML</sub>", "M<sub>ML</sub>"),
#     Type = "LRT"
#   ) |>
#   dplyr::select(Model, Type, Estimator, pval)

# Histogram alternative to the p-value distribution

fig2 <- null_dat |>
  filter(iter %in% 1:20) |>
  filter(convergence == 0) |>
  dplyr::select(model, iter, betacat, delta, target, term, p_wald, p_satt) |>
  pivot_longer(cols = p_wald:p_satt, names_to = "type", values_to = "pval") |>
  mutate(
    Estimator = if_else(grepl("_reml", model), "REML", "ML"),
    model = gsub("_reml", "", model),
    Model = if_else(model == "beta", "&beta;", "M"),
    Model = paste0(Model, "<sub>", Estimator, "</sub>"),
    Type = if_else(type == "p_satt", "Satterthwaite DDF", "Wald p-values"),
  ) |>
  filter(!is.na(pval)) |>
  dplyr::select(Model, Type, Estimator, pval) |>

  mutate(
    Type = factor(Type, levels = c("Wald p-values", "Satterthwaite DDF"))
  ) |>

  ggplot(aes(pval)) +
  geom_histogram(aes(y = after_stat(density)), boundary = 0) +

  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +

  facet_wrap(Model ~ Type, ncol = 5) +
  theme_classic() +

  labs(x = "*p*-value", y = "Density") +

  theme(
    axis.text = element_text(size = axis_text_size),
    strip.text = ggtext::element_markdown(
      hjust = 0,
      margin = margin(t = 1, b = 1, unit = "pt")
    ),
    axis.text.x = element_text(hjust = c(0, 0.5, 1)),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = ggtext::element_markdown(size = axis_text_size)
  )


# Extraction of error rates are aggregated over several steps

s1 <- null_dat |>
  filter(convergence == 0) |>
  dplyr::select(model, iter, betacat, delta, target, term, p_wald, p_satt) |>
  pivot_longer(cols = p_wald:p_satt, names_to = "type", values_to = "pval") |>

  filter(!is.na(pval)) |>

  dplyr::summarise(
    .by = c(model, betacat, iter, term, type),
    n = dplyr::n(),
    sig = sum(pval < 0.05),
    error = sig / n
  )


# Collapse strata and terms within an iteration: the permutation is the
# resampling unit, so one error rate per iteration per model.
s1_iter <- s1 |>
  dplyr::summarise(.by = c(model, iter, type), error = mean(error))

model_error_type1 <- s1_iter |>
  dplyr::summarise(
    .by = c(model, type),
    m.error = mean(error),
    s.error = stats::sd(error),
    n_iter = dplyr::n(),
    mcse = s.error / sqrt(n_iter)
  ) |>
  dplyr::mutate(lo = m.error - 2.576 * mcse, hi = m.error + 2.576 * mcse)

t1e <- s1_iter |>
  dplyr::summarise(
    .by = c(model, type),
    m.error = mean(error),
    s.error = stats::sd(error),
    n_iter = dplyr::n(),
    mcse = s.error / sqrt(n_iter)
  )


model_error_type1_betacat <- s1 |>
  dplyr::summarise(.by = c(model, iter, type, betacat), error = mean(error)) |>
  dplyr::summarise(
    .by = c(model, betacat, type),
    m.error = mean(error),
    s.error = stats::sd(error),
    n_iter = dplyr::n(),
    mcse = s.error / sqrt(n_iter)
  ) |>
  dplyr::mutate(lo = m.error - 2.576 * mcse, hi = m.error + 2.576 * mcse)


# Only colour strips in x-direction
strip <- ggh4x::strip_themed(
  background_x = lapply(
    colors,
    \(x) element_rect(fill = scales::alpha(x, 0.3), colour = NA)
  )
)

fig3 <- model_error_type1_betacat |>

  mutate(
    betacat = case_when(
      betacat == "xlow" ~ "0 - 0.2",
      betacat == "low" ~ "0.2 - 0.4",
      betacat == "mid" ~ "0.4 - 0.6",
      betacat == "high" ~ "0.6 - 0.8",
      betacat == "xhigh" ~ "0.8 - 1"
    ),
    model = dplyr::recode(
      model,
      beta = "&beta;<sub>ML</sub>",
      beta_reml = "&beta;<sub>REML</sub>",
      m = "M<sub>ML</sub>",
      m_reml = "M<sub>REML</sub>"
    ),
    type = dplyr::recode(
      type,
      p_wald = "Wald p-values",
      p_satt = "Satterthwaite DDF"
    )
  ) |>

  ggplot(aes(model, 100 * m.error, shape = type, fill = type)) +

  geom_errorbar(
    aes(ymin = 100 * lo, ymax = 100 * hi),
    width = 0.5,
    position = position_dodge(width = 0.3)
  ) +
  geom_hline(yintercept = 5, lty = 2, color = "gray") +

  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  ggh4x::facet_wrap2(~betacat, ncol = 5, strip = strip) +
  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.text.x = ggtext::element_markdown(size = axis_text_size),
    legend.position = "right",
    legend.position.inside = c(0.75, 0.15),
    legend.background = element_rect(),
    legend.text = element_text(size = 8),
    legend.direction = "vertical",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = axis_title_size)
  ) +

  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    breaks = c(2.5, 5, 7.5, 10),
    limits = c(4, 10),
    expand = c(0, 0)
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(y = "Type I error rate (%)")


# 03. Standard errors per category and inferential procedure ##################

se_dat <- perm_dat |>
  dplyr::filter(convergence == 0, term %in% alt_terms, delta == 0) |>
  dplyr::mutate(
    divisor = dplyr::if_else(grepl("^beta", model), log(2), 1),
    se_M = std.error / divisor,
    # Satterthwaite can return a degenerate denominator degrees of
    # freedom on a fit whose participant variance sits on the zero
    # boundary (@fig-re-sd-failedsatterthwaite). The inner replacement
    # keeps qt() off those values rather than letting it warn.
    ddf_ok = !is.na(ddf) & ddf > 0,
    # The multiplier that each test applies to the standard error.
    p_wald = stats::qnorm(0.975),
    p_satt = dplyr::if_else(
      ddf_ok,
      stats::qt(0.975, dplyr::if_else(ddf_ok, ddf, 1)),
      NA_real_
    )
  ) |>
  dplyr::select(model, iter, target, betacat, se_M, ddf, p_wald, p_satt) |>
  tidyr::pivot_longer(
    cols = c(p_wald, p_satt),
    names_to = "type",
    values_to = "crit"
  ) |>
  dplyr::mutate(
    type = factor(type, levels = c("p_wald", "p_satt")),
    # The half-width of the rejection region in M units.
    halfwidth = crit * se_M
  )

# The permutation is the resampling unit, as it is everywhere else: a site
# level summary is formed within an iteration first and the iterations are
# then averaged. The within-iteration summary is a median because the
# standard error of the beta arm is right skewed in the extreme strata, where
# a few sites with near-boundary fitted means dominate a mean.
se_stratum <- se_dat |>
  dplyr::filter(type == "p_wald") |>
  dplyr::summarise(
    .by = c(model, betacat, iter),
    med_se = stats::median(se_M, na.rm = TRUE)
  )


# Annotation
se_annotation <- se_stratum |>
  dplyr::mutate(model = lab_model(model), betacat = lab_betacat(betacat)) |>
  filter(iter == 1, betacat == "0 - 0.2") |>
  mutate(med_se = if_else(grepl("REML", model), med_se + 0.013, med_se - 0.005))


fig4 <- se_stratum |>
  dplyr::mutate(model = lab_model(model), betacat = lab_betacat(betacat)) |>

  ggplot(aes(model, med_se, fill = model)) +

  geom_violin() +

  # geom_point(
  #   size = 1.5,
  #   shape = 21,
  #   fill = "white",
  #   position = position_jitter(width = 0.2)
  # ) +
  ggh4x::facet_wrap2(~betacat, ncol = 5, strip = strip) +
  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title.x = ggtext::element_markdown(),
    legend.position = "none",
    strip.text = element_blank(),

    axis.text.x = element_blank(),
    axis.title.y = ggtext::element_markdown(size = axis_title_size),
    axis.ticks.x.bottom = element_blank()
  ) +

  ggtext::geom_richtext(
    data = se_annotation,
    aes(model, med_se, label = model),
    fill = NA,
    label.color = NA, # remove background and outline
    label.padding = grid::unit(rep(0, 4), "pt"), # remove padding
    size = 2.5
  ) +

  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0.09, 0.15), expand = 0) +

  scale_x_discrete(sec.axis = dup_axis()) +
  labs(x = NULL, y = "Standard error<br>(M-value units)")


# 04. Power (only in the REML-Satterthwaite procedure) ########################

# The nominal level.
alpha <- 0.05

# Monte Carlo critical value.
#
# The functions first sorts the p-values and then
# sets then uses the alpha * length(p) as the critical value.
crit_value <- function(p, alpha) {
  p <- sort(p[!is.na(p)])
  k <- floor(alpha * length(p))
  # Too few draws to place a rejection region at this level: reject nothing.
  if (k < 1) return(0)
  p[k]
}

# One row per (model, test type, iteration, site). We use the wald-type
# contrasts and the data from the LRT test (ML only).
#
# `is_alt` is ground truth. These are the sites with added effect.
power_dat <- perm_dat |>
  dplyr::filter(convergence == 0, term %in% alt_terms) |>
  dplyr::select(
    model,
    iter,
    target,
    betacat,
    delta,
    delta_abs,
    estimate,
    p_wald,
    p_satt
  ) |>
  tidyr::pivot_longer(
    cols = c(p_wald, p_satt),
    names_to = "type",
    values_to = "pval"
  ) |>

  # Applied to nulls and alternatives together: dropping the Satterthwaite
  # failures from one side only would change which sites a cell is scored
  # on. The retained counts are reported per cell in @tbl-power-calibration.
  dplyr::filter(!is.na(pval)) |>
  dplyr::mutate(
    is_alt = delta != 0,
    sgn = sign(delta),
    model = factor(model),
    type = factor(type, levels = c("p_wald", "p_satt"))
  )

stopifnot(any(power_dat$is_alt), any(!power_dat$is_alt))


# Calibration, pooled over strata.
crit_global <- power_dat |>
  dplyr::filter(!is_alt) |>
  dplyr::summarise(
    .by = c(model, type),
    n_null = dplyr::n(),
    # Actual size at the nominal cut-off.
    size_nominal = mean(pval <= alpha),
    c_hat = crit_value(pval, alpha),
    # By construction this should sit at alpha. It is kept
    # as a check on the rejection region, not as a result.
    size_adjusted = mean(pval <= c_hat)
  )

# Calibration within stratum.
crit_strata <- power_dat |>
  dplyr::filter(!is_alt) |>
  dplyr::summarise(
    .by = c(model, type, betacat),
    n_null = dplyr::n(),
    size_nominal = mean(pval <= alpha),
    c_hat = crit_value(pval, alpha),
    size_adjusted = mean(pval <= c_hat)
  )

# The rejection region must be conservative in every cell.
stopifnot(
  all(crit_global$size_adjusted <= alpha + 1e-6),
  all(crit_strata$size_adjusted <= alpha + 1e-6)
)

# Per-iteration power in the finest cell of the design: model, test type,
# stratum, effect magnitude and direction. Three calibrations are carried side
# by side so that the effect of the adjustment is readable.
power_cell <- power_dat |>
  dplyr::filter(is_alt) |>
  dplyr::inner_join(
    crit_global |> dplyr::select(model, type, c_global = c_hat),
    by = c("model", "type")
  ) |>
  dplyr::inner_join(
    crit_strata |> dplyr::select(model, type, betacat, c_stratum = c_hat),
    by = c("model", "type", "betacat")
  ) |>
  dplyr::summarise(
    .by = c(model, type, betacat, delta_abs, sgn, iter),
    n_sites = dplyr::n(),
    pow_nominal = mean(pval <= alpha),
    pow_adj_global = mean(pval <= c_global),
    pow_adj_stratum = mean(pval <= c_stratum)
  )

# Aggregate power_cell to a requested grouping. The mean is taken over the
# finer cells within an iteration first, so that an unbalanced cell (caused by
# a convergence failure) does not weight the pooled estimate.
power_agg <- function(dat, by = character()) {
  dat |>
    dplyr::summarise(
      .by = dplyr::all_of(c("model", "type", "iter", by)),
      dplyr::across(dplyr::starts_with("pow_"), mean)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("pow_"),
      names_to = "calibration",
      values_to = "power"
    ) |>
    dplyr::summarise(
      .by = dplyr::all_of(c("model", "type", "calibration", by)),
      m.power = mean(power),
      s.power = stats::sd(power),
      n_iter = dplyr::n(),
      mcse = s.power / sqrt(n_iter)
    ) |>
    dplyr::mutate(lo = m.power - 1.96 * mcse, hi = m.power + 1.96 * mcse)
}

power_overall <- power_agg(power_cell)
power_delta <- power_agg(power_cell, by = c("delta_abs"))
power_strata <- power_agg(power_cell, by = c("betacat", "delta_abs"))
power_sign <- power_agg(power_cell, by = c("delta_abs", "betacat", "sgn"))


fig5 <- power_sign |>
  dplyr::filter(calibration == "pow_adj_stratum") |>
  dplyr::mutate(
    Model = lab_model(model),
    type = lab_type(type),
    betacat = lab_betacat(betacat)
  ) |>

  filter(type == "Satterthwaite DDF", grepl("_reml", model)) |>

  ggplot(aes(
    delta_abs * sgn,
    100 * m.power,
    color = betacat,
    group = paste(Model, betacat)
  )) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2, alpha = 0.5) +
  ggh4x::facet_grid2(~Model) +
  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = ggtext::element_markdown(size = axis_title_size),
    axis.title.x = ggtext::element_markdown(size = axis_title_size),
    axis.title.y = ggtext::element_markdown(size = axis_title_size),
    axis.text.x = element_text(size = axis_text_size),
    strip.text = ggtext::element_markdown(hjust = 0),
    legend.key.spacing.y = unit(-2, 'mm'),
    strip.background = element_blank(),
  ) +
  scale_colour_manual(values = colors) +
  scale_shape_manual(values = c(21, 24, 23, 22)) +
  scale_x_continuous(
    breaks = c(-0.8, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.8),
    labels = c("-.8", "-.4", "", "", "", "", ".4", ".8")
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 5)) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Added effect (&delta;, M-value)",
    y = "Size-adjusted<br>power (%)"
  )


# 05. Results from full run #################################################

# The results directory is VERSIONED in the same way as the permutation
# directory: full_v2/ holds the four-arm run (beta, beta_reml, m, m_reml)
# with `p_wald` and `p_satt` per contrast, the latter on the REML arms only.
# The earlier two-arm files in the derived_data root (beta-model-full.RDS,
# m-model-full.RDS) carry a different summary schema and are superseded by it.
full_dir <- here::here("analysis/data/derived_data/full_v3")
full_arms <- c("beta_reml", "m_reml")
full_files <- file.path(full_dir, sprintf("%s.RDS", full_arms))

# Unlike the permutation study, a missing file is NOT refitted on the fly: the
# four arms together are close to a day of wall-clock time (see the timings
# reported below), so what is missing is reported instead.
if (!all(file.exists(full_files)))
  stop(
    "missing arm files in full_v2/: ",
    paste(basename(full_files[!file.exists(full_files)]), collapse = ", "),
    ". Run analysis/R/methylation-case-full.R first.",
    call. = FALSE
  )

# Optionally trim the read to a subset of fixed effects. All terms are kept by
# default; the four time contrasts plus the intercept over 770 441 sites is
# ~3.9 million rows per arm, so setting this to the contrast(s) actually
# reported is the lever if the render is memory-bound.
full_keep_terms <- NULL

# One arm at a time. Each arm object is on the order of 10^2 MB, so the
# per-fit diagnostics are joined on WITHIN the arm and the arm is dropped
# before the next one is read, rather than binding four whole arms and joining
# afterwards.
read_full_arm <- function(file) {
  arm <- readRDS(file)

  summaries <- if (is.null(full_keep_terms)) arm$summaries else
    dplyr::filter(arm$summaries, term %in% full_keep_terms)

  # Carry convergence, singularity and the participant SD across from the
  # evaluation slot, exactly as perm_dat does above: the arms fail at
  # different rates, so a failed fit has to be droppable from the tests
  # rather than quietly changing which sites each arm is scored on.
  dat <- summaries |>
    dplyr::inner_join(
      dplyr::select(arm$evaluations, target, convergence, singular, re_sd),
      by = "target"
    ) |>
    # Constant within an arm, so stored as factors: at this number of
    # rows the same three labels as character vectors are hundreds of
    # MB of pointers for no extra information.
    dplyr::mutate(
      model = factor(arm$model, levels = full_arms),
      estimator = factor(
        if (isTRUE(arm$meta$reml)) "REML" else "ML",
        levels = c("ML", "REML")
      ),
      scale = factor(arm$meta$scale, levels = c("beta", "M"))
    )

  # seqwrap_summarise() drops failed targets rather than returning empty
  # rows, so a target present in one slot and not the other would silently
  # shrink the arm here.
  stopifnot(nrow(dat) == nrow(summaries))

  evaluations <- dplyr::mutate(
    arm$evaluations,
    model = factor(arm$model, levels = full_arms)
  )

  # `class` and `message` come back as character vectors named by target;
  # the names are dropped so the bound frame does not carry four sets of
  # them.
  errors <- if (NROW(arm$errors)) {
    arm$errors |>
      dplyr::mutate(
        class = unname(class),
        message = unname(message),
        model = factor(arm$model, levels = full_arms)
      )
  } else NULL

  out <- list(
    dat = dat,
    evaluations = evaluations,
    errors = errors,
    meta = arm$meta
  )

  rm(arm, summaries)
  gc()

  out
}

full_list <- lapply(full_files, read_full_arm)

# full_dat:    one row per (model, site, fixed effect) -- the frame the
#              per-position analysis is read from.
# full_evals:  one row per (model, site) -- fit-level diagnostics.
# full_errors: one row per condition raised (warnings and errors), by stage.
# full_meta:   one row per arm -- what was fitted, on how many cores, how long.
full_dat <- dplyr::bind_rows(lapply(full_list, `[[`, "dat"))
full_evals <- dplyr::bind_rows(lapply(full_list, `[[`, "evaluations"))
full_errors <- dplyr::bind_rows(lapply(full_list, `[[`, "errors"))
full_meta <- dplyr::bind_rows(lapply(full_list, `[[`, "meta"))

rm(full_list)
gc()

# `reml` in the evaluation slot is read off the fit object, not off the arm
# name (see eval_fun in the script), so this is the check that the arm
# labelling agrees with what was actually estimated. The arms are also all
# fitted on the same set of sites, which is what makes them comparable
# site-by-site.
stopifnot(
  setequal(as.character(full_evals$model), full_arms),
  all(full_evals$reml == grepl("_reml$", as.character(full_evals$model))),
  length(unique(full_meta$k)) == 1L,
  length(unique(full_meta$n)) == 1L
)

# The number of sites and samples the full run was fitted on.
n_full_sites <- unique(full_meta$k)
n_full_samples <- unique(full_meta$n)

# The contrasts: 1_baseline is the reference level, so the time terms are the
# differences from baseline. The intercept is present in full_dat but is not a
# contrast and is excluded here.
full_terms <- setdiff(unique(full_dat$term), "(Intercept)")


# Combine the full results for the target models
res <- full_dat |>
  filter(estimator == "REML") |>
  filter(
    term %in%
      c("time2_acute", "time3_loading", "time4_unloading", "time5_reloading")
  ) |>
  mutate(est_M = if_else(grepl("beta", model), estimate / log(2), estimate))

d <- res

density_fun <- function(d, TERM) {
  dd <- filter(d, term == TERM) |>
    pull(est_M) |>
    density(from = -1, to = 1)

  data.frame(x = dd$x, y = dd$y) |>
    mutate(time = TERM)
}


results1 <- bind_rows(
  density_fun(res, "time2_acute"),
  density_fun(res, "time3_loading"),
  density_fun(res, "time4_unloading"),
  density_fun(res, "time5_reloading")
) |>

  mutate(
    time = case_when(
      time == "time2_acute" ~ "Acute",
      time == "time3_loading" ~ "Loading",
      time == "time4_unloading" ~ "Unloading",
      time == "time5_reloading" ~ "Reloading"
    ),
    time = factor(
      time,
      levels = c("Acute", "Loading", "Unloading", "Reloading")
    )
  ) |>

  ggplot(aes(x, y)) +
  geom_line() +
  facet_wrap(~time, ncol = 4) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-.8", "-.4", 0, ".4", ".8")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = axis_text_size),
    strip.background = element_blank(),
    axis.title.x = ggtext::element_markdown(size = axis_title_size),
    strip.text = element_text(hjust = 0, size = axis_title_size)
  ) +
  labs(x = "Estimated effects (M-value units)")


# full_errors |>
#   filter(stage == "fit") |>
#   summarise(.by = c(type, model), n = n()) |>
#   mutate(percentage = n / n_full_sites) |>
#   print()
#
# full_errors |>
#   distinct(message)
#
#
# res |>
#
#   mutate(.by = c(model, term), fdr = p.adjust(p_satt, method = "fdr")) |>
#
#   filter(p_satt < 0.05) |>
#
#   mutate(sign = if_else(estimate > 0, "up", "down")) |>
#   summarise(.by = c(term, model, sign), n = n()) |>
#   mutate(identity = if_else(sign == "up", n, -n)) |>
#   ggplot(aes(term, identity, fill = model)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   geom_hline(yintercept = 0)
#

# 06. Timing experiment #####################################################

# The annotation file
gset <- readRDS(
  here::here("analysis/data/derived_data/seaborne-gset-normalized.RDS")
)
anno <- minfi::getAnnotation(gset)
nsites <- nrow(anno)

rm(gset)
gc()

# The timimg results
timing_dir <- here::here("analysis/data/derived_data/timing_v1")

timing_cells <- readRDS(file.path(timing_dir, "timing-cells.RDS"))
timing_scaling <- readRDS(file.path(timing_dir, "timing-scaling.RDS"))

# Size block: is the cost linear in k?
#
# `timing-scaling.RDS` carries an `exponent` and a `linear` flag from a log-log
# fit of elapsed time on k. That fit has NO intercept term, so it tests whether
# elapsed time is PROPORTIONAL to k, which cannot hold when there is a fixed
# per-call cost of about a minute: at the k used here the intercept dominates,
# and the exponent it returns (~0.16 to 0.35) measures that dominance rather
# than the scaling. The saved `linear` column is therefore FALSE for all four
# arms and is not read here.
#
# The test with the intercept restored is elapsed = a + b * k^c, with linearity
# at c = 1, and the curvature check is whether a quadratic term is needed. Both
# are computed from the raw cells.
timing_size <- dplyr::filter(timing_cells, block == "size")

timing_fit <- timing_size |>
  dplyr::group_by(arm) |>
  dplyr::group_modify(function(d, key) {
    lin <- stats::lm(elapsed ~ k, data = d)
    quad <- stats::lm(elapsed ~ k + I(k^2), data = d)

    nl <- tryCatch(
      stats::nls(
        elapsed ~ a + b * k^c,
        data = d,
        start = list(
          a = stats::coef(lin)[[1]],
          b = stats::coef(lin)[[2]],
          c = 1
        ),
        control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
      ),
      error = function(e) NULL
    )

    # Profile intervals on an nls fit can fail to converge; the Wald
    # interval is the fallback, and NA is reported rather than a number
    # of unknown provenance if both fail.
    ci <- if (is.null(nl)) c(NA_real_, NA_real_) else {
      s <- summary(nl)$coefficients["c", ]
      prof <- tryCatch(
        suppressMessages(stats::confint(nl, "c")),
        error = function(e) NULL,
        warning = function(w) NULL
      )
      if (is.null(prof)) {
        s[["Estimate"]] + c(-1.96, 1.96) * s[["Std. Error"]]
      } else {
        as.numeric(prof)
      }
    }

    data.frame(
      r2 = summary(lin)$r.squared,
      exponent = if (is.null(nl)) NA_real_ else stats::coef(nl)[["c"]],
      exponent_lo = ci[1],
      exponent_hi = ci[2],
      p_quadratic = summary(quad)$coefficients["I(k^2)", 4]
    )
  }) |>
  dplyr::ungroup()


timing_tab <- timing_scaling$scaling |>
  dplyr::select(arm, scale, estimator, overhead_s, per_target_s) |>
  dplyr::inner_join(timing_fit, by = "arm") |>
  dplyr::inner_join(
    dplyr::select(
      timing_scaling$extrapolation,
      arm,
      k_full,
      predicted_h,
      observed_h,
      error_pct
    ),
    by = "arm"
  ) |>
  dplyr::mutate(per_target_ms = per_target_s * 1000)

# Approximate start up cost
# This calculates the intercept for each arm and averages

start_up <- timing_cells |>
  dplyr::filter(block == "size") |>

  group_by(arm) |>
  summarise(intercept = coef(lm(elapsed ~ k))[1]) |>
  ungroup() |>
  summarise(intercept = mean(intercept)) |>
  pull(intercept)

# Time figure 1, a regression line shows the cost with increasing number of
# m (number of targets)
time1 <- timing_cells |>
  dplyr::filter(block == "size") |>

  mutate(
    arm = case_when(
      arm == "beta" ~ "&beta;<sub>ML</sub>",
      arm == "m" ~ "M<sub>ML</sub>",
      arm == "beta_reml" ~ "&beta;<sub>REML</sub>",
      arm == "m_reml" ~ "M<sub>REML</sub>"
    ),
    arm = forcats::fct_reorder(arm, elapsed),
    arm = forcats::fct_rev(arm)
  ) |>

  ggplot(aes(k, elapsed, fill = arm, color = arm)) +
  geom_point(alpha = 1, size = 3, color = "black", shape = 21) +

  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +

  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  labs(x = "m (number of targets)", y = "Elapsed time\n(seconds)") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 300), expand = c(0, 0)) +

  geom_hline(yintercept = start_up, lty = 2, color = "gray") +
  annotate(
    "text",
    x = 3000,
    y = start_up + 15,
    label = "Start-up cost",
    color = "gray40",
    size = 2.5
  ) +

  theme(
    axis.text.x = element_text(size = axis_text_size),
    legend.text = ggtext::element_markdown(),
    legend.position = "none",
    axis.title = element_text(size = axis_title_size),
  )


# Obs and pred for top costing model
annotation_time <- timing_tab |>
  filter(arm == "beta_reml") |>
  dplyr::select(observed_h, predicted_h) |>
  pivot_longer(cols = everything()) |>
  mutate(
    name = gsub("_h", "", name),
    name = stringr::str_to_sentence(name),
    value = if_else(value == max(value), value + 0.1, value - 0.3)
  )

# Timimng (and prediction) of the full experiment
time2 <- timing_tab |>
  dplyr::mutate(model = lab_model(arm)) |>
  dplyr::select(
    model,
    estimator,
    overhead_s,
    per_target_ms,
    r2,
    exponent,
    exponent_lo,
    exponent_hi,
    p_quadratic,
    predicted_h,
    observed_h,
    error_pct
  ) |>

  dplyr::select(model, estimator, predicted_h, observed_h) |>
  pivot_longer(cols = predicted_h:observed_h) |>
  mutate(
    sites = "m = 770441",
    name = if_else(name == "predicted_h", "Predicted", "Observed"),
    Model = if_else(name == "Predicted", model, "")
  ) |>
  mutate(
    model = forcats::fct_reorder(model, value),
    model = forcats::fct_rev(model)
  ) |>

  ggplot(aes(x = sites, value, fill = model, shape = name)) +
  geom_point(size = 3) +

  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    limits = c(2, 12),
    expand = c(0, 0),
    sec.axis = dup_axis()
  ) +

  labs(x = "", y = "Elapsed time\n(hours)") +

  ggtext::geom_richtext(
    aes(label = Model, color = NULL),
    position = position_nudge(x = 0.2),
    fill = NA,
    label.color = NA, # remove background and outline
    label.padding = grid::unit(rep(0, 4), "pt"), # remove padding)
    size = 2.5,
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = axis_text_size),
    legend.position = "none",
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.line.y.left = element_blank(),
    axis.title.y.left = element_blank(),
    axis.text.y.right = element_text(size = axis_text_size),
    axis.title.y.right = element_text(size = axis_text_size)
  ) +
  annotate(
    "text",
    x = 0.65,
    y = annotation_time$value + c(0.5, -0.5),
    label = annotation_time$name,
    size = 2.5,
    color = "gray40"
  ) +
  annotate(
    "curve",
    y = annotation_time$value[1] + 0.3,
    x = 0.8,
    yend = annotation_time$value[1] + c(-0.1),
    xend = 0.95,
    curvature = c(-0.2),
    arrow = arrow(type = "closed", length = unit(2, "mm")),
    color = "gray40"
  ) +

  annotate(
    "curve",
    y = annotation_time$value[2] - 0.3,
    x = 0.8,
    yend = annotation_time$value[2] + c(0.1),
    xend = 0.95,
    curvature = c(0.2),
    arrow = arrow(type = "closed", length = unit(2, "mm")),
    color = "gray40"
  )


# 07. Combine full figure  ##################################################

figure5 <- cowplot::plot_grid(
  cowplot::plot_grid(
    fig1,
    fig3,
    fig4,
    align = "v",
    axis = "lr",
    ncol = 1,
    rel_heights = c(0.4, 1, 0.6, 0.8)
  ),
  cowplot::plot_grid(
    NULL,
    fig5,
    results1,
    rel_widths = c(0.1, 1, 1),
    ncol = 3,
    hjust = 0.8
  ),
  cowplot::plot_grid(
    NULL,
    time1,
    time2,
    NULL,
    align = "h",
    axis = "lr",
    rel_widths = c(0.1, 0.9, 0.6, 0.1),
    ncol = 4
  ),
  ncol = 1,
  rel_heights = c(1.5, 0.5, 0.5)
) +
  annotate(
    "text",
    x = c(0.03, 0.03, 0.03, 0.03, 0.5, 0.03, 0.6),
    y = c(0.98, 0.85, 0.57, 0.37, 0.37, 0.2, 0.2),
    label = c("A", "B", "C", "D", "E", "F", "G")
  )


saveRDS(figure5, here::here("analysis/figures/figure-5.RDS"))


ggsave(
  here::here("analysis/figures/figure-5.pdf"),
  figure5,
  device = cairo_pdf,
  height = 180,
  width = 170,
  units = "mm"
)

ggsave(
  here::here("analysis/figures/figure-5.png"),
  figure5,
  bg = "white",
  dpi = 600,
  height = 180,
  width = 170,
  units = "mm"
)
