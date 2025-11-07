
# to summarise brms models

# sum_fun_brms <- function(x){
#
#   cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$fixed))),
#                              coef(summary(x))$fixed,
#
#                              row.names = NULL)
#
#   return(cond_effects)
#
# }


# Summary function for brms::brm models with posterior probabilities and Bayes factors
sum_fun_brms <- function(x) {
  # Extract fixed effects summary
  fixef_mat <- brms::fixef(x, summary = TRUE)
  cond_effects <- data.frame(
    term = rownames(fixef_mat),
    fixef_mat,
    row.names = NULL
  )

  # Compute posterior probabilities (P(β > 0))
  post_draws <- as.data.frame(brms::posterior_samples(x, pars = rownames(fixef_mat)))
  post_probs <- sapply(post_draws, function(p) mean(p > 0))
  cond_effects$posterior_prob_gt0 <- post_probs[cond_effects$term]

  # Compute Bayes factors (against null = 0)
  bf_list <- brms::bayes_factor(x, hypothesis = paste0(rownames(fixef_mat), " = 0"))
  cond_effects$bayes_factor <- bf_list$bf

  return(cond_effects)
}



sum_fun_brms1 <- function(x) {
  fe <- fixef(x)
  cond_effects <- data.frame(
    coef = rownames(fe),
    fe,
    row.names = NULL
  )
  return(cond_effects)
}


sum_fun_lmer <- function(x){

  # Get name from x
  # geneid <- names(x)

  cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))))),
                             coef(summary(x)),

                             row.names = NULL) #%>%
  #  mutate(geneid = geneid)

  return(cond_effects)
}

sum_fun_lmer1 <- function(model) {
  tryCatch({
    coefs <- summary(model)$coefficients
    data.frame(
      term = rownames(coefs),
      coefs,
      row.names = NULL
    )
  }, error = function(e) {
    message("Error in sum_fun_lmer: ", e$message)
    return(NULL)
  })
}


# summary function for lme
sum_fun_lme <- function(x) {
  # Extract the fixed effects table
  ttab <- summary(x)$tTable

  # Create a data frame with rownames as a column
  cond_effects <- data.frame(
    term = rownames(ttab),
    ttab,
    row.names = NULL
  )

  return(cond_effects)
}




# summary function for nlme::lme
sum_fun_lme1 <- function(
    x,
    conf.int = FALSE,
    level = 0.95,
    ci.method = c("t", "intervals"),
    drop_intercept = FALSE
) {
  if (!inherits(x, "lme")) {
    stop("This function expects an 'nlme::lme' model object.")
  }

  # Extract and validate tTable
  tt <- summary(x)$tTable
  if (is.null(tt)) stop("Could not find tTable in summary(x).")

  # Build base data.frame with term column
  df <- data.frame(
    term = rownames(tt),
    as.data.frame(tt),
    row.names = NULL,
    check.names = FALSE
  )

  # Standardize column names if present
  col_map <- c(
    "Value" = "estimate",
    "Std.Error" = "std.error",
    "DF" = "df",
    "t-value" = "statistic",
    "p-value" = "p.value"
  )
  intersecting <- intersect(names(col_map), names(df))
  names(df)[match(intersecting, names(df))] <- unname(col_map[intersecting])

  # Optional: drop intercept row
  if (isTRUE(drop_intercept)) {
    df <- df[df$term != "(Intercept)", , drop = FALSE]
  }

  # Confidence intervals
  if (isTRUE(conf.int)) {
    ci.method <- match.arg(ci.method)
    if (!all(c("estimate", "std.error") %in% names(df))) {
      stop("Needed columns not found to compute confidence intervals.")
    }

    if (ci.method == "t") {
      # Classic t-based Wald CI using per-term df if available
      alpha <- 1 - level
      if ("df" %in% names(df) && !anyNA(df$df)) {
        crit <- stats::qt(1 - alpha/2, df = df$df)
      } else {
        # Fallback: use large-sample normal approx
        crit <- stats::qnorm(1 - alpha/2)
      }
      df$conf.low  <- df$estimate - crit * df$std.error
      df$conf.high <- df$estimate + crit * df$std.error

    } else if (ci.method == "intervals") {
      # Use nlme::intervals for fixed effects
      ints <- try(nlme::intervals(x, which = "fixed", level = level), silent = TRUE)
      if (inherits(ints, "try-error")) {
        warning("nlme::intervals() failed; falling back to t-based intervals.")
        alpha <- 1 - level
        if ("df" %in% names(df) && !anyNA(df$df)) {
          crit <- stats::qt(1 - alpha/2, df = df$df)
        } else {
          crit <- stats::qnorm(1 - alpha/2)
        }
        df$conf.low  <- df$estimate - crit * df$std.error
        df$conf.high <- df$estimate + crit * df$std.error
      } else {
        fx <- as.data.frame(ints$fixed)
        # columns are typically: lower, est., upper
        fx$term <- rownames(ints$fixed)
        rownames(fx) <- NULL
        names(fx) <- sub("^est\\.$", "estimate_int", names(fx)) # avoid clash
        # Join by term while preserving original order
        m <- match(df$term, fx$term)
        df$conf.low  <- fx$lower[m]
        df$conf.high <- fx$upper[m]
        # (We keep df$estimate from tTable to match testing columns)
      }
    }
  }

  # Reorder columns nicely if present
  desired <- c("term", "estimate", "std.error", "df", "statistic", "p.value",
               "conf.low", "conf.high")
  df <- df[, intersect(desired, names(df)), drop = FALSE]

  df
}

#evaluation function for lme

evalfun_lme <- function(x) {
  # A function that returns a data frame with the p-values from
  # the residual diagnostics
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  results <- data.frame(
    pval.unif = unif$p-value,
  )

  return(results)
}


# summary and evaluation functions for GLMM models
sumfun <- function(x) {
  # A function that returns a data frame with the coefficients
  # from the fitted model
  cond_effects <- data.frame(
    cbind(
      data.frame(
        coef = rownames(coef(summary(x))$cond)
      )
    ),
    coef(summary(x))$cond,
    row.names = NULL
  )

  return(cond_effects)
}

# Updated sumfun
sumfun2 <- function(x) {
  cond_summary <- coef(summary(x))$cond

  # Convert to data frame and add rownames as a column
  cond_effects <- as.data.frame(cond_summary)
  cond_effects$coef <- rownames(cond_summary)

  # Reorder columns so 'coef' comes first
  cond_effects <- cond_effects[, c("coef", setdiff(names(cond_effects), "coef"))]

  return(cond_effects)
}



evalfun <- function(x) {

  sim <- DHARMa::simulateResiduals(x, n = 1000)

  disp <- DHARMa::testDispersion(sim, plot = FALSE)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  zinfl <- DHARMa::testZeroInflation(sim, plot = FALSE)

  results <- data.frame(pval.disp = disp$p.value,
                        pval.unif = unif$p.value,
                        pval.zinfl = zinfl$p.value)

  return(results)
}

evalfun_brms <- function(x) {
  # Posterior predictive checks
  pp <- brms::pp_check(x, type = "stat", plot = FALSE)

  # Extract basic fit statistics
  loo_res <- brms::loo(x)
  waic_res <- brms::waic(x)

  results <- data.frame(
    loo_elpd = loo_res$estimates["elpd_loo", "Estimate"],
    loo_se = loo_res$estimates["elpd_loo", "SE"],
    waic = waic_res$estimates["waic", "Estimate"],
    waic_se = waic_res$estimates["waic", "SE"]
  )

  return(results)
}
