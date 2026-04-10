################################################################################
# Prepare observed data for simulations
#
# 01. Packages
# IF m1 result do not exixts ---
# 02. Download and prep data
# 03. A preliminary model
# ------------------------------
# 04. Load data
#
# Notes: Results from model m1 is used to generate simulated data, i.e. it must
# be present to sun the simulation script. The script checks if m1 is available,
# otherwise it re-runs the estimations process (time-consuming).
#
################################################################################

# 01. Packages ################################################################

source("R/check-packages.R")

library(seqwrap)
library(tidyverse)
library(ggtext)
library(glmmTMB)
library(edgeR)
library(DESeq2)


# 02. Download and prep data #################################################

if (!dir.exists("data/")) dir.create("data/")

if (!file.exists("data/m1_results.RDS")) {
  dat <- seqwrap::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  # Calculate observed library sizes
  observed_mean_libsize <- mean(libsize$lib.size) / 20
  observed_cv_libsize <- sd(libsize$lib.size) / mean(libsize$lib.size)

  # 03. A preliminary model #####################################################

  # A preliminary model is fitted using a conditional NB distribution.
  # The purpose of the conditional model is to estimate distributions
  # of parameter estimates. We will use these for simulations

  # A summary function to return the dispersion parameter with SE on the log scale
  # mean(predict(x, type = "link)) will give us the predicted log counts.
  # We will put this in the eval fun to also get estimates of the parameters in
  # the generic summary function.
  sigma_summary <- function(x) {
    out <- data.frame(
      dispersion = data.frame(summary(x$sdr))["betadisp", 1],
      dispersion.se = data.frame(summary(x$sdr))["betadisp", 2],
      log_mu = mean(predict(x, type = "link"))
    )
    return(out)
  }

  m1 <- seqwrap_compose(
    data = countdat,
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = glmmTMB,
    eval_fun = sigma_summary,
    arguments = list(
      formula = y ~ time * group + offset(ln_efflibsize) + (1 | id),
      family = glmmTMB::nbinom2
    )
  )

  m1_results <- seqwrap(
    m1,
    return_models = FALSE,

    cores = parallel::detectCores()
  )
  saveRDS(m1_results, "data/m1_results.RDS")
}


# 04. Load data  ##############################################################

m1_results <- readRDS("data/m1_results.RDS")


m1_sum <- seqwrap_summarise(m1_results, verbose = FALSE)

# Fit a trend to the dispersion data from m1
# save the data in a convenient format.
dispersion_dat <- m1_sum$evaluations

# Fit a loess model, this model can be used in the
# simulation function to represent the mean-dispersion
# relationship
trend_model_observed <- loess(
  dispersion ~ log_mu,
  data = dispersion_dat,
  span = 0.7,
  weights = 1 / (dispersion.se^2)
)

trend_model_observed_noweights <- loess(
  dispersion ~ log_mu,
  data = dispersion_dat,
  span = 0.7
)


new_data <- data.frame(log_mu = seq(from = 0, to = 14, by = 0.1))


# Trend model predictions
dispersion_pred <- data.frame(
  log_mu = seq(from = 0, to = 14, by = 0.1),
  pred = predict(trend_model_observed, newdata = new_data)
) |>
  mutate(sd = trend_model_observed$s)


# # Check the dispersion fit
dispersion_dat |>
  ggplot(aes(log_mu, dispersion)) +
  geom_point(alpha = 0.2) +

  geom_ribbon(
    data = dispersion_pred,
    aes(log_mu, pred, ymin = pred - sd, ymax = pred + sd),
    alpha = 0.2
  ) +

  geom_line(data = dispersion_pred, aes(log_mu, pred), color = "red")

# # Overall distributions of observed effects
# m1_sum$summaries |>
#   ggplot(aes(estimate)) +
#   geom_density() +
#   facet_wrap(~ term, scales = "free")
#
# # Random effect distribution on the log-normal scale
# m1_sum$summaries |>
#   filter(term == "sd__(Intercept)") |>
#   summarise(m = mean(log(estimate)),
#             s = sd(log(estimate)))
#
