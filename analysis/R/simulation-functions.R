###############################################################################
#
# Paper-specific simulation helpers
#
# These functions are sourced by analysis/paper/make-docs.R to enable
# running the simulation study.
# Functions depend on:
# tidyverse, edgeR, glmmTMB, DESeq2, lmerTest, DHARMa,
# broom.mixed, and the seqwrap packages.
#
###############################################################################

#' Summary function for dispersion parameter
#'
#' A summary function to return the dispersion parameter with SE on the log
#' scale mean(predict(x, type = "link)) will give us the predicted log counts.
#' We will put this in the eval fun to also get estimates of the parameters in
#' the generic summary function.
#'
#' @param x A fitted glmmTMB model object.
sigma_summary <- function(x) {
  if (is.null(x$fit$convergence)) {
    conv <- 1
  } else {
    conv <- x$fit$convergence[[1]]
  }

  out <- data.frame(
    dispersion = data.frame(summary(x$sdr))["betadisp", 1],
    dispersion.se = data.frame(summary(x$sdr))["betadisp", 2],
    log_mu = mean(stats::predict(x, type = "link")),
    convergence = conv,
    pdHess = x$sdr$pdHess
  )
  return(out)
}


#' A summary function for the lmer model of transformed counts
#' it will return a the singularity diagnostics from lme4.
#'
#' @param x A fitted lmer model object.
lmer_summary <- function(x) {
  out <- data.frame(isSingular = lme4::isSingular(x))

  return(out)
}


#' A summary function for the Poisson models. This will use the convergence
#' diagnostic in glmmTMB to indicate convergence. The pdHess indicator from the
#' Hessian matrices
#' see https://stackoverflow.com/questions/79110546/glmmtmb-convergence-messages
#'
#' @param x A fitted glmmTMB Poisson model object.
poisson_summary <- function(x) {
  if (is.null(x$fit$convergence)) {
    conv <- 1
  } else {
    conv <- x$fit$convergence[[1]]
  }

  out <- data.frame(convergence = conv, pdHess = x$sdr$pdHess)
  return(out)
}


#' Sigma summary common for both Poisson and NB models on real data
#'
#' A summary function used in modelling
#'
#' @param x A fitted glmmTMB model object.
sigma_summary2 <- function(x) {
  if (is.null(x$fit$convergence)) {
    conv <- 1
  } else {
    conv <- x$fit$convergence[[1]]
  }

  ### Simulating scaled resid ###
  # This creates simulated residuals for test of uniformity and dispersion.
  # See
  # https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
  # for details.
  residObj <- DHARMa::simulateResiduals(x, plot = FALSE)

  # Saving tests from DHARMa
  unif <- DHARMa::testUniformity(residObj, plot = FALSE)
  disp <- DHARMa::testDispersion(residObj, plot = FALSE)

  # Saving SD of Obs Level Rand Eff in the Poisson model

  tidy_coef <- broom.mixed::tidy(x)

  if (any(tidy_coef$group %in% "seq_sample_id")) {
    olre.sd <- tidy_coef |>
      dplyr::filter(group == "seq_sample_id") |>
      dplyr::pull(estimate)
  } else {
    olre.sd <- NA
  }

  # Combine all in a data frame
  out <- data.frame(
    dispersion = data.frame(summary(x$sdr))["betadisp", 1],
    dispersion.se = data.frame(summary(x$sdr))["betadisp", 2],
    olre.sd = olre.sd,
    log_mu = mean(stats::predict(x, type = "link")),
    convergence = conv,
    pdHess = x$sdr$pdHess,
    aic = stats::AIC(x),
    unif.p = unif$p.value,
    unif.stat = unif$statistic,
    disp.p = disp$p.value,
    disp.stat = disp$statistic
  )
  return(out)
}


#' Download simulations and methylation case-study results
#'
#' The simulations and the methylation case study in this project are
#' computationally intensive. Their results are available at
#' Chidimma Echebiri; Ellefsen, Stian; Ahmad, Rafi; Hammarstrom, Daniel,
#' 2026, "Simulated data sets for: seqwrap: an R package for flexible iterative
#'  fitting of high-dimensional data", https://doi.org/10.18710/I7U71O,
#'  DataverseNO.
#'  This function downloads all files in the data set from dataverse.
#'
#' The data set holds two kinds of files, told apart by their directory label
#' on Dataverse:
#'
#' * Simulation results (`estimates/`, `estimates2/`, `evaluations/`,
#'   `evaluations2/`, `simdata/`, `simdata2/`) and the data set README go to
#'   `dest_dir`, i.e. `analysis/data/raw_data/`.
#' * Methylation case-study results carry the prefix `derived_data/`
#'   (`derived_data/permutation_v4/`, `derived_data/full_v3/`,
#'   `derived_data/full_v2/`, `derived_data/timing_v1/` and the two array
#'   files in `derived_data/` itself). The prefix is stripped and the files go
#'   to `derived_dir`, i.e. `analysis/data/derived_data/`, which is where the
#'   scripts in `analysis/R/` would have written them.
#'
#' @param doi The DOI in our data set.
#' @param dest_dir The destination folder for simulation results. Should be
#' placed in raw data (added to .gitignore)
#' @param derived_dir The destination folder for files whose directory label
#' starts with `derived_data/`.
#' @param server For reuse of the function, if other dataverse servers
#' are to be used.
#' @param overwrite If we need to overwrite individual files.
download_dataverse <- function(
  doi = "doi:10.18710/I7U71O",
  dest_dir = here::here("analysis/data/raw_data"),
  derived_dir = here::here("analysis/data/derived_data"),
  server = "dataverse.no",
  overwrite = FALSE
) {
  # Clean DOI format
  doi <- sub("https://doi.org/", "doi:", doi, fixed = TRUE)
  doi <- sub("http://doi.org/", "doi:", doi, fixed = TRUE)
  if (!grepl("^doi:", doi)) doi <- paste0("doi:", doi)

  # Get file metadata
  meta_url <- sprintf(
    "https://%s/api/datasets/:persistentId/?persistentId=%s",
    server,
    doi
  )
  meta_file <- tempfile(fileext = ".json")

  system2(
    "curl.exe",
    args = c("-k", "-L", shQuote(meta_url), "-o", shQuote(meta_file))
  )

  meta <- jsonlite::fromJSON(meta_file, simplifyDataFrame = FALSE)
  files <- meta$data$latestVersion$files
  file.remove(meta_file)

  message("Found ", length(files), " files to download")

  purrr::walk(files, function(f) {
    dir <- if (is.null(f$directoryLabel)) "" else f$directoryLabel

    # Route by directory label: the methylation case-study files are labelled
    # derived_data/... on Dataverse and belong under analysis/data/derived_data
    # (the prefix is stripped); everything else is a simulation file and goes
    # to dest_dir as before.
    if (grepl("^derived_data(/|$)", dir)) {
      root <- derived_dir
      dir <- sub("^derived_data/?", "", dir)
    } else {
      root <- dest_dir
    }
    out_dir <- if (nzchar(dir)) file.path(root, dir) else root
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    dest <- file.path(out_dir, f$dataFile$filename)

    if (file.exists(dest) && !overwrite) {
      message(
        "Skipping (already exists): ",
        file.path(dir, f$dataFile$filename)
      )
      return(invisible(NULL))
    }

    file_url <- sprintf(
      "https://%s/api/access/datafile/%s",
      server,
      f$dataFile$id
    )
    message("Downloading: ", file.path(dir, f$dataFile$filename))
    system2(
      "curl.exe",
      args = c("-k", "-L", shQuote(file_url), "-o", shQuote(dest))
    )
  })

  invisible(dest_dir)
}


#' Extract simulation results
#'
#' This function extract simulations results for analysis. Simulation results
#' are either downloaded (download_dataverse) or simulated using the simulation
#' wrapper.
#'
#' @param evaluations_path Path to evaluation result files.
#' @param estimates_path Path to estimate result files.
#' @param populationeffects_path Path to population-effect result files.
#' @param disp_scenario Marker for the simulation scenario (dispersion setting).
extract_simulations <- function(
  evaluations_path = here::here("analysis/data/raw_data/evaluations"),
  estimates_path = here::here("analysis/data/raw_data/estimates"),
  populationeffects_path = here::here("analysis/data/raw_data/popeffect"),
  disp_scenario = "s1"
) {
  # Evaluations

  # The evaluations can be combined despite having different shape. Non-
  # available columns will be NA.

  eval_files <- list.files(evaluations_path)
  evaluations <- list()
  for (i in seq_along(eval_files)) {
    evaluations[[i]] <- readRDS(file.path(evaluations_path, eval_files[i])) |>
      dplyr::mutate(file = eval_files[i], disp_scenario = disp_scenario)
  }

  evaluations <- dplyr::bind_rows(evaluations)

  # Estimates
  est_files <- list.files(estimates_path)
  estimates <- list()
  for (i in seq_along(est_files)) {
    estimates[[i]] <- readRDS(file.path(estimates_path, est_files[i])) |>
      dplyr::mutate(file = est_files[i], disp_scenario = disp_scenario)
  }

  estimates <- dplyr::bind_rows(estimates)

  # Population effects
  pop_files <- list.files(populationeffects_path)
  popeffects <- list()
  for (i in seq_along(pop_files)) {
    popeffects[[i]] <- readRDS(file.path(
      populationeffects_path,
      pop_files[i]
    )) |>
      dplyr::mutate(file = pop_files[i], disp_scenario = disp_scenario)
  }

  popeffects <- dplyr::bind_rows(popeffects)

  return(list(
    populationeffects = popeffects,
    estimates = estimates,
    evaluations = evaluations
  ))
}


#' Filter function
#'
#' A function for filter by expression, some genes will have low expression
#' due to sampling variability, these are removed in this function.
#' The function is used in the simulation wrappers.
#'
#' counts, a count data frame. First column are target id.
#' metadata, a metadat data frame containing group/time combinations
filter_fun <- function(counts, metadata) {
  ## Filter by expression
  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    counts[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(metadata$condition, metadata$time)
  )

  counts_filtered <- counts[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(counts_filtered[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    tibble::rownames_to_column(var = "seq_sample_id") |>
    dplyr::select(-group)

  # Combine all meta data
  metadata <- metadata |>
    dplyr::inner_join(libsize, by = "seq_sample_id") |>
    dplyr::mutate(
      efflibsize = (lib.size * norm.factors) /
        stats::median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  return(list(counts = counts_filtered, metadata = metadata))
}


#' A (wrapper) function for simulating data sets
#'
#' @param nullgenes the number of genes with null effects
#' @param condB_true number of genes with non-zero effects in baseline group diffs
#' @param condB_timet2_true number of genes with non-zero effects in interaction
#' @param dispersion_model The model used for dispersion estimates
#' @param dataset data set id.
simulate_datasets <- function(
  nullgenes = 7500,
  condB_true = 1250,
  condB_time2_true = 1250,
  dispersion_model = NULL,
  dataset
) {
  ngenes <- nullgenes + condB_true + condB_time2_true

  ## Set fixed (population level) effects
  beta0 <- stats::runif(
    nullgenes + condB_true + condB_time2_true,
    min = 1.5,
    max = 7
  )
  conditionB <- c(
    rep(0, nullgenes),
    stats::runif(
      condB_true,
      min = 0.2,
      max = 1
    ) *
      sample(c(-1, 1), condB_true, prob = c(0.5, 0.5), replace = TRUE),
    rep(0, condB_time2_true)
  )

  timet2 <- stats::rnorm(nullgenes + condB_true + condB_time2_true, 0, 0.1)
  timet3 <- stats::rnorm(nullgenes + condB_true + condB_time2_true, 0, 0.2)

  conditionB_timet2 <- rep(
    0,
    nullgenes +
      condB_true +
      condB_time2_true
  )

  conditionB_timet3 <- c(
    rep(0, nullgenes),
    rep(0, condB_true),
    stats::runif(
      condB_true,
      min = 0.2,
      max = 1
    ) *
      sample(c(-1, 1), condB_true, prob = c(0.25, 0.75), replace = TRUE)
  )

  # Simulate random effects
  # approximately based on observed data
  b0_values <- stats::rlnorm(
    nullgenes + condB_true + condB_time2_true,
    meanlog = -2.07,
    sdlog = 1
  )

  # Set b1 and b2 distribution to ~small
  b1_values <- rep(0, nullgenes + condB_true + condB_time2_true)
  b2_values <- rep(0, nullgenes + condB_true + condB_time2_true)

  # Simulate data #
  simdat <- seqwrap::simcounts2(
    n1 = 32,
    n2 = 40,
    beta0 = beta0,
    conditionB = conditionB,
    timet2 = timet2,
    timet3 = timet3,
    conditionB_timet2 = conditionB_timet2,
    conditionB_timet3 = conditionB_timet3,
    b0 = b0_values,
    b1 = b1_values,
    b2 = b2_values,
    # Using the trend model from observed data
    phi_model = dispersion_model,
    lib_size_mean = 10^6,
    lib_size_cv = 0.145,
    max_prop = 0.02
  )

  # Subdivide data sets into different sample sizes
  # small = 8 + 10
  # medium = 16 + 20 (similar to the observed)
  # large = 32 + 40

  metadata_small <- simdat$metadata |>
    dplyr::filter(id %in% c(paste0("A", 1:8), paste0("B", 1:10)))
  metadata_medium <- simdat$metadata |>
    dplyr::filter(id %in% c(paste0("A", 1:16), paste0("B", 1:20)))
  metadata_large <- simdat$metadata |>
    dplyr::filter(id %in% c(paste0("A", 1:32), paste0("B", 1:40)))

  counts_small <- simdat$counts |>
    dplyr::select(gene, tidyselect::all_of(metadata_small$seq_sample_id))

  counts_medium <- simdat$counts |>
    dplyr::select(gene, tidyselect::all_of(metadata_medium$seq_sample_id))

  counts_large <- simdat$counts |>
    dplyr::select(gene, tidyselect::all_of(metadata_large$seq_sample_id))

  # Filter low expression genes, the resulting count tables contain
  # all genes that are to be used in simulations.

  combined_data <- list(
    small = filter_fun(counts = counts_small, metadata = metadata_small),
    medium = filter_fun(counts = counts_medium, metadata = metadata_medium),
    large = filter_fun(counts = counts_large, metadata = metadata_large)
  )

  # Genes present in data sets after filtering
  genes_small <- combined_data[[1]]$counts$gene
  genes_medium <- combined_data[[2]]$counts$gene
  genes_large <- combined_data[[3]]$counts$gene

  filtered_genes <- data.frame(
    size = c(
      rep("small", length(genes_small)),
      rep("medium", length(genes_medium)),
      rep("large", length(genes_large))
    ),
    target = c(genes_small, genes_medium, genes_large)
  ) |>
    tidyr::expand_grid(term = c("conditionB", "timet3:conditionB"))

  # Save population effects
  population_effects <- dplyr::bind_rows(
    data.frame(
      target = 1:ngenes,
      term = rep("conditionB", ngenes),
      population_effect = conditionB,
      dataset = dataset
    ),
    data.frame(
      target = 1:ngenes,
      term = rep("timet3:conditionB", ngenes),
      population_effect = conditionB_timet3,
      dataset = dataset
    )
  ) |>
    tidyr::expand_grid(size = c("small", "medium", "large")) |>
    dplyr::inner_join(filtered_genes)

  return(list(
    simdat = simdat,
    combined_data = combined_data,
    population_effects = population_effects
  ))
}


#' Model 1 and 2 simulation function
#'
#' Model 1 and 2 are the naive and informed negative binomial models
#' weighted_loess, should a weighted loess be used for mean-dispersion estimates
#' dataset, when used in a loop dataset indicate the index
#' dofit, if false only data management is done
m1_m2_sim <- function(
  combined_data,
  dataset,
  dofit = TRUE,
  weighted_loess = TRUE,
  CORES = parallel::detectCores()
) {
  evaluations <- list()
  summaries <- list()
  evaluations2 <- list()
  summaries2 <- list()

  for (k1 in seq_along(combined_data)) {
    ms1 <- seqwrap::seqwrap_compose(
      data = combined_data[[k1]]$counts, # These are the filtered counts
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = sigma_summary,
      targetdata = NULL,
      arguments = list(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1 | id),
        family = glmmTMB::nbinom2
      )
    )

    if (dofit) {
      ms1_results <- seqwrap::seqwrap(
        ms1,
        return_models = FALSE,
        # subset = 1:1100,
        verbose = FALSE,
        cores = CORES
      )

      ms1_sum <- seqwrap::seqwrap_summarise(ms1_results, verbose = FALSE)

      evaluations[[k1]] <- ms1_sum$evaluations |>
        dplyr::mutate(
          model = "m1",
          datasets = dataset,
          size = names(combined_data)[k1]
        )

      summaries[[k1]] <- ms1_sum$summaries |>
        dplyr::mutate(
          model = "m1",
          datasets = dataset,
          size = names(combined_data)[k1]
        )
    }

    ## Model 2 ##
    # get successful targets
    targets <- evaluations[[k1]] |>
      dplyr::filter(!is.na(dispersion.se)) |>
      dplyr::distinct(target) |>
      dplyr::pull(target)

    ## Fitting a model for the mean-dispersion relationship

    # Fit a trend to the dispersion data from m1
    # save the data in a convenient format. Using the log_mu_raw (average
    # observed counts) allows for a prior on dispersion also for genes with
    # unsuccessful fits in model 1. First we gather all dispersion data.

    dispersion_dat <- evaluations[[k1]] |>
      dplyr::filter(target %in% targets)

    # Calculate the raw observed counts from the data
    raw_log_counts <- data.frame(
      target = as.character(combined_data[[k1]]$counts[, 1]),
      log_mu_raw = log(rowMeans(combined_data[[k1]]$counts[, -1]))
    )

    # Adding log raw counts to the dispersion df for modeling.
    dispersion_dat <- dispersion_dat |>
      dplyr::inner_join(raw_log_counts)

    # Fit a loess model, using log_mu_raw as the predictor
    if (weighted_loess) {
      trend_model <- stats::loess(
        dispersion ~ log_mu_raw,
        data = dispersion_dat,
        span = 0.7,
        weights = 1 / (dispersion.se^2)
      )
    } else {
      trend_model <- stats::loess(
        dispersion ~ log_mu_raw,
        data = dispersion_dat,
        span = 0.7
      )
    }

    # Predict dispersion for each gene based on log raw counts
    # and combine into a prior.

    dispersion_prior <- data.frame(
      gene = combined_data[[k1]]$counts[, 1],
      pred = round(
        stats::predict(
          trend_model,
          newdata = data.frame(
            log_mu_raw = log(
              rowMeans(
                combined_data[[k1]]$counts[, -1]
              )
            )
          )
        ),
        3
      ),
      s = round(
        trend_model$s,
        3
      )
    ) |>
      dplyr::mutate(prior = paste0("normal(", pred, ",", s, ")"))

    # Gather all distributions of estimates that can be used as prior information
    # in subsequent model. Fixed effects are centered on 0. We will use the SD
    # for priors

    estimate_distributions <- summaries[[k1]] |>
      dplyr::filter(target %in% targets) |>
      dplyr::select(target, term, estimate) |>

      dplyr::summarise(
        .by = term,
        m = mean(estimate),
        s = stats::sd(estimate)
      ) |>
      dplyr::filter(!term %in% c("(Intercept)", "sd__(Intercept)"))

    # Extract the random effects distribution to fit a gamma distribution
    random_sd_estimate <- summaries[[k1]] |>
      dplyr::filter(target %in% targets, term == "sd__(Intercept)") |>
      dplyr::select(target, term, estimate) |>
      dplyr::pull(estimate)

    # The gamma distribution is parameterized using a shape and a rate
    # parameter. It looks like this prior will lead to a push towards 0,
    # consider adding a constant to push away from zero...
    # TODO this may needs testing.
    mean_sd <- mean(random_sd_estimate)
    var_sd <- stats::var(random_sd_estimate)
    shape_param <- 2 # mean_sd^2 / var_sd

    # Here we prepare priors for the fixed effects, in this version all fixed
    # effects, except the intercept, will have regularizing priors corresponding to
    # the distributions of effects seen in the naive models.
    Priors_df <- dplyr::bind_rows(
      data.frame(
        prior = paste0("normal(0,", round(estimate_distributions$s, 2), ")"),
        class = rep("fixef", 5),
        coef = estimate_distributions$term
      ),
      data.frame(
        prior = paste0(
          "gamma(",
          round(mean_sd, 2),
          ",",
          2,
          ")"
        ),
        class = "ranef",
        coef = "id"
      )
    )

    # We want to use the mean-dispersion relationship to add a prior for the
    # dispersion parameter. This means that we need a gene specific prior

    Priors_list <- list()
    for (j in 1:nrow(combined_data[[k1]]$counts)) {
      df <- dplyr::bind_rows(
        Priors_df,
        data.frame(
          prior = dispersion_prior[
            dispersion_prior$gene == combined_data[[k1]]$counts[j, 1],
            4
          ],
          class = "fixef_disp",
          coef = "1"
        )
      )

      Priors_list[[j]] <- df
    }

    # seqwrap accepts target-wise data frames as a list,
    # this makes it easier to specify target-specific priors.

    # Here we specify priors based the results from ms1

    ms2 <- seqwrap::seqwrap_compose(
      data = combined_data[[k1]]$counts,
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = sigma_summary,
      targetdata = Priors_list,
      arguments = alist(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1 | id),
        family = glmmTMB::nbinom2,
        priors = data.frame(
          prior = prior,
          class = class,
          coef = coef
        )
      )
    )

    ms2_results <- seqwrap::seqwrap(
      ms2,
      return_models = FALSE,
      verbose = FALSE,
      #  subset = 1:50,
      cores = CORES
    )

    ms2_sum <- seqwrap::seqwrap_summarise(ms2_results, verbose = FALSE)

    evaluations2[[k1]] <- ms2_sum$evaluations |>
      dplyr::mutate(
        model = "m2",
        datasets = dataset,
        size = names(combined_data)[k1]
      )

    summaries2[[k1]] <- ms2_sum$summaries |>
      dplyr::mutate(
        model = "m2",
        datasets = dataset,
        size = names(combined_data)[k1]
      )
  }

  return(list(
    evaluations_m1 = dplyr::bind_rows(evaluations),
    evaluations_m2 = dplyr::bind_rows(evaluations2),
    summaries_m1 = dplyr::bind_rows(summaries),
    summaries_m2 = dplyr::bind_rows(summaries2)
  ))
}


#' Model 3 simulation function
#'
#' This is the model for transformed counts data the function performs
#' transformation and fits models over a combined_data object.
m3_sim <- function(combined_data, dataset, dofit = TRUE,
                   CORES = parallel::detectCores()) {
  evaluations <- list()
  summaries <- list()

  for (k3 in seq_along(combined_data)) {
    # The VST transformation

    # Safe check of any NA's in the data
    safe_counts <- combined_data[[k3]]$counts[
      stats::complete.cases(combined_data[[k3]]$counts[, -1]),
      -1
    ]

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = safe_counts,
      colData = combined_data[[k3]]$metadata,
      design = ~ time * condition
    )

    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    vst_mat <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(
      dds,
      blind = FALSE,
      fitType = "parametric"
    ))
    vst_dat <- cbind(data.frame(
      gene = combined_data[[k3]]$counts[
        stats::complete.cases(combined_data[[k3]]$counts[, -1]),
        1
      ],
      as.data.frame(vst_mat)
    ))

    ms3 <- seqwrap::seqwrap_compose(
      data = vst_dat,
      metadata = combined_data[[k3]]$metadata,
      samplename = "seq_sample_id",
      modelfun = lmerTest::lmer,
      eval_fun = lmer_summary,
      arguments = list(
        formula = y ~ time * condition + (1 | id)
      )
    )

    if (dofit) {
      ms3_results <- seqwrap::seqwrap(
        ms3,
        return_models = FALSE,
        verbose = FALSE,
        #   subset = 1:50,
        cores = CORES
      )

      ms3_sum <- seqwrap::seqwrap_summarise(ms3_results, verbose = FALSE)

      evaluations[[k3]] <- ms3_sum$evaluations |>
        dplyr::mutate(
          model = "m3",
          datasets = dataset,
          size = names(combined_data)[k3]
        )
      summaries[[k3]] <- ms3_sum$summaries |>
        dplyr::mutate(
          model = "m3",
          datasets = dataset,
          size = names(combined_data)[k3]
        )
    }
  }

  evaluations <- dplyr::bind_rows(evaluations)
  summaries <- dplyr::bind_rows(summaries)

  return(list(
    evaluations_m3 = dplyr::bind_rows(evaluations),
    summaries_m3 = dplyr::bind_rows(summaries)
  ))
}


#' Model 4 and 5 simulation function
#'
#' This function does the Poisson models with observation level random effects
#' both in a naive and informed version.
m4_m5_sim <- function(combined_data, dataset, dofit = TRUE,
                      CORES = parallel::detectCores()) {
  evaluations <- list()
  summaries <- list()
  evaluations2 <- list()
  summaries2 <- list()

  for (k1 in seq_along(combined_data)) {
    ms1 <- seqwrap::seqwrap_compose(
      data = combined_data[[k1]]$counts, # These are the filtered counts
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = poisson_summary,
      targetdata = NULL,
      arguments = list(
        formula = y ~
          time *
            condition +
            offset(ln_efflibsize) +
            (1 | id) +
            (1 | seq_sample_id),
        family = stats::poisson
      )
    )

    if (dofit) {
      ms1_results <- seqwrap::seqwrap(
        ms1,
        return_models = FALSE,
        verbose = FALSE,
        # subset = 1:1100,
        cores = CORES
      )

      ms1_sum <- seqwrap::seqwrap_summarise(ms1_results, verbose = FALSE)

      evaluations[[k1]] <- ms1_sum$evaluations |>
        dplyr::mutate(
          model = "m4",
          datasets = dataset,
          size = names(combined_data)[k1]
        )

      summaries[[k1]] <- ms1_sum$summaries |>
        dplyr::mutate(
          model = "m4",
          datasets = dataset,
          size = names(combined_data)[k1]
        )
    }

    ## Model 2 ##
    # get successful targets
    targets <- evaluations[[k1]] |>
      dplyr::filter(convergence == 0) |>
      dplyr::distinct(target) |>
      dplyr::pull(target)

    estimate_distributions <- summaries[[k1]] |>
      dplyr::filter(target %in% targets) |>
      dplyr::select(target, term, group, estimate) |>

      dplyr::summarise(
        .by = c(term, group),
        m = mean(estimate),
        s = stats::sd(estimate)
      ) |>
      dplyr::filter(!term %in% c("(Intercept)", "sd__(Intercept)"))

    # Extract the random effects distribution to fit a gamma distribution
    # on the participant level intercepts
    random_sd_estimate <- summaries[[k1]] |>
      dplyr::filter(
        target %in% targets,
        term == "sd__(Intercept)",
        group == "id"
      ) |>
      dplyr::select(target, term, estimate) |>
      dplyr::pull(estimate)

    random_sd_estimate_obs <- summaries[[k1]] |>
      dplyr::filter(
        target %in% targets,
        term == "sd__(Intercept)",
        group == "seq_sample_id"
      ) |>
      dplyr::select(target, term, estimate) |>
      dplyr::pull(estimate)

    # The gamma distribution is parameterized using a shape and a rate
    # parameter. It looks like this prior will lead to a push towards 0,
    # consider adding a constant to push away from zero...
    # TODO this may needs testing.
    mean_sd <- mean(random_sd_estimate)
    var_sd <- stats::var(random_sd_estimate)
    shape_param <- 2 # mean_sd^2 / var_sd
    # Observation level
    mean_sd_obs <- mean(random_sd_estimate_obs)
    var_sd_obs <- stats::var(random_sd_estimate_obs)
    shape_param_obs <- 2 # mean_sd^2 / var_sd

    # Here we prepare priors for the fixed effects, in this version all fixed
    # effects, except the intercept, will have regularizing priors corresponding to
    # the distributions of effects seen in the naive models.

    # The priors for the random effects are added, glmmTMB accepts coefficient
    # index.
    Priors_df <- dplyr::bind_rows(
      data.frame(
        prior = paste0("normal(0,", round(estimate_distributions$s, 2), ")"),
        class = rep("fixef", 5),
        coef = estimate_distributions$term
      ),
      data.frame(
        prior = paste0(
          "gamma(",
          c(round(mean_sd, 2), round(mean_sd_obs, 2)),
          ",",
          2,
          ")"
        ),
        class = "ranef",
        coef = c("1", "2")
      )
    )

    # We want to use the mean-dispersion relationship to add a prior for the
    # dispersion parameter. This means that we need a gene specific prior

    Priors_list <- list()
    for (j in 1:nrow(combined_data[[k1]]$counts)) {
      Priors_list[[j]] <- Priors_df
    }

    # seqwrap accepts target-wise data frames as a list,
    # this makes it easier to specify target-specific priors.

    # Here we specify priors based the results from ms1

    ms2 <- seqwrap::seqwrap_compose(
      data = combined_data[[k1]]$counts,
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = poisson_summary,
      targetdata = Priors_list,
      arguments = alist(
        formula = y ~
          time *
            condition +
            offset(ln_efflibsize) +
            (1 | id) +
            (1 | seq_sample_id),
        family = stats::poisson,
        priors = data.frame(
          prior = prior,
          class = class,
          coef = coef
        )
      )
    )

    ms2_results <- seqwrap::seqwrap(
      ms2,
      return_models = FALSE,
      verbose = FALSE,
      #  subset = 1:50,
      cores = CORES
    )

    ms2_sum <- seqwrap::seqwrap_summarise(ms2_results, verbose = FALSE)

    evaluations2[[k1]] <- ms2_sum$evaluations |>
      dplyr::mutate(
        model = "m5",
        datasets = dataset,
        size = names(combined_data)[k1]
      )

    summaries2[[k1]] <- ms2_sum$summaries |>
      dplyr::mutate(
        model = "m5",
        datasets = dataset,
        size = names(combined_data)[k1]
      )
  }

  return(list(
    evaluations_m4 = dplyr::bind_rows(evaluations),
    evaluations_m5 = dplyr::bind_rows(evaluations2),
    summaries_m4 = dplyr::bind_rows(summaries),
    summaries_m5 = dplyr::bind_rows(summaries2)
  ))
}


#' Extract raw counts in simulated data
#'
#' This function calculates mean counts for each data set in the simulation
#' at the relevant terms.
#'
#' @param sim_folder Path to the simulation folder
extract_rawcounts <- function(
  sim_folder = here::here("analysis/data/raw_data/simdata2/raw/")
) {
  # Sample ids for each data set
  samps <- list(
    small_baseline_samp = c(
      paste0(paste0("A", 1:8), "_t1"),
      paste0(paste0("B", 1:10), "_t1")
    ),
    small_time3_samp = c(
      paste0(paste0("A", 1:8), "_t3"),
      paste0(paste0("B", 1:10), "_t3")
    ),
    medium_baseline_samp = c(
      paste0(paste0("A", 1:16), "_t1"),
      paste0(paste0("B", 1:20), "_t1")
    ),
    medium_time3_samp = c(
      paste0(paste0("A", 1:16), "_t3"),
      paste0(paste0("B", 1:20), "_t3")
    ),
    large_baseline_samp = c(
      paste0(paste0("A", 1:32), "_t1"),
      paste0(paste0("B", 1:40), "_t1")
    ),
    large_time3_samp = c(
      paste0(paste0("A", 1:32), "_t3"),
      paste0(paste0("B", 1:40), "_t3")
    )
  )

  # A small meta data for the data sets
  meta_data <- data.frame(
    size = c("small", "small", "medium", "medium", "large", "large"),
    term = rep(c("conditionB", "timet3:conditionB"), 3)
  )

  # Extract the observed mean/variances
  mean_var <- function(samps, counts) {
    temp <- counts |>
      dplyr::select(gene, dplyr::matches(samps))

    out <- data.frame(
      target = counts[, 1],
      m = rowMeans(temp[, -1]),
      var = apply(temp[, -1], 1, stats::var)
    )

    return(out)
  }

  files <- list.files(sim_folder)

  d <- list()
  for (i in seq_along(files)) {
    dataset <- gsub(".RDS", "", gsub("dataset_", "", files[i]))
    temp <- readRDS(paste0(sim_folder, files[i]))

    dsub <- list()
    for (j in 1:nrow(meta_data)) {
      dsub[[j]] <- mean_var(samps[[j]], temp$counts) |>
        dplyr::mutate(
          size = meta_data[j, 1][[1]],
          term = meta_data[j, 2][[1]],
          dataset = dataset
        )
    }

    d[[i]] <- dplyr::bind_rows(dsub)
  }

  out <- dplyr::bind_rows(d)
  return(out)
}


#' Cheap PSOCK worker start-up for a block of seqwrap() calls
#'
#' `seqwrap()` creates a PSOCK cluster and destroys it again on every call, and
#' each worker runs the project `.Rprofile` -- so renv activates once per
#' worker, per call. Measured on a 12-core machine that is roughly 19 s of
#' fixed cost per `seqwrap()` call. `sim_wrap1()` and `sim_wrap2()` each make
#' 15 calls per dataset (3 sample sizes x 5 models) across 10 datasets, i.e.
#' 150 calls, so a large part of each wrapper is spent starting workers rather
#' than fitting models.
#'
#' Pointing `R_PROFILE_USER` at an empty file makes the workers skip
#' `.Rprofile`, and pushing `.libPaths()` through `R_LIBS` keeps them on the
#' same renv library so package versions are unchanged. This takes the fixed
#' cost to about 4 s per call.
#'
#' Because the workers no longer run `.Rprofile`, renv does not hide the base R
#' library from them: they get the renv library and sandbox first (identical to
#' this session) with the system library appended as a fallback. A package that
#' is MISSING from the renv library would therefore be taken silently from the
#' system library on a worker. `fast_workers_begin()` asserts that the packages
#' doing the work resolve to the same installation in parent and worker.
#'
#' `normalizePath()` matters here: for a package that is *attached* in this
#' session `find.package()` reports the path it was loaded from, and renv's
#' library entries are junctions into the renv cache -- so the parent would say
#' ".../renv/cache/.../seqwrap" where a worker, resolving through
#' `.libPaths()`, says ".../renv/library/.../seqwrap". Same installation, two
#' spellings. Comparing the resolved directory removes that false alarm while
#' still catching a genuine fallback to the system library.
#'
#' Pair every call with `fast_workers_end()`, ideally through `on.exit()` so
#' that the session is restored even if the caller aborts.
#'
#' @param pkgs Packages that must resolve identically on the workers. These are
#'   the ones referenced by the model functions and by the summary and
#'   evaluation functions used in the simulations.
#' @param quiet Suppress the pbapply progress bar, which `seqwrap()` prints on
#'   every call and which is noise inside a loop.
#' @return A state object to hand back to `fast_workers_end()`.
fast_workers_begin <- function(
  pkgs = c("seqwrap", "glmmTMB", "broom.mixed", "lmerTest", "lme4"),
  quiet = TRUE
) {
  state <- list(
    env = Sys.getenv(c("R_PROFILE_USER", "R_LIBS"), unset = NA),
    profile = tempfile(fileext = ".R"),
    pbo = NULL
  )
  file.create(state$profile)

  Sys.setenv(
    R_PROFILE_USER = state$profile,
    R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)
  )
  if (quiet) state$pbo <- pbapply::pboptions(type = "none")

  stamp <- function(pk) {
    vapply(
      pk,
      function(p) {
        paste(
          normalizePath(find.package(p), winslash = "/", mustWork = FALSE),
          utils::packageVersion(p)
        )
      },
      character(1)
    )
  }

  parent <- stamp(pkgs)
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  worker <- parallel::clusterCall(cl, stamp, pkgs)

  if (!all(vapply(worker, identical, logical(1), parent))) {
    fast_workers_end(state)
    stop(
      "Worker package resolution differs from the parent session. ",
      "A package is probably missing from the renv library and would be ",
      "taken from the system library on the workers.",
      call. = FALSE
    )
  }

  state
}


#' Restore the session after `fast_workers_begin()`
#'
#' @param state The value returned by `fast_workers_begin()`.
fast_workers_end <- function(state) {
  if (!is.null(state$pbo)) pbapply::pboptions(state$pbo)

  for (v in names(state$env)) {
    if (is.na(state$env[[v]])) {
      Sys.unsetenv(v)
    } else {
      do.call(Sys.setenv, stats::setNames(list(state$env[[v]]), v))
    }
  }
  unlink(state$profile)
  invisible(NULL)
}


#' Simulation wrapper 1
#'
#' A wrapper for the simulation functions for scenario 1
#'
#' @param cores Number of cores to be used
#' @param seed Set the seed.
#' @param overwrite Should available results be overwritten?
sim_wrap1 <- function(cores, seed = 1, overwrite = FALSE) {
  set.seed(1)

  # 150 seqwrap() calls follow (10 datasets x 3 sample sizes x 5 models);
  # see fast_workers_begin() for why worker start-up dominates without this.
  .fw <- fast_workers_begin()
  on.exit(fast_workers_end(.fw), add = TRUE)

  raw_dir <- here::here("analysis/data/raw_data/simdata/raw")
  clean_dir <- here::here("analysis/data/raw_data/simdata/clean")
  popeff_dir <- here::here("analysis/data/raw_data/simdata/popeffect")
  est_dir <- here::here("analysis/data/raw_data/estimates")
  eval_dir <- here::here("analysis/data/raw_data/evaluations")

  for (i in 1:10) {
    d <- simulate_datasets(
      nullgenes = 7500,
      condB_true = 1250,
      condB_time2_true = 1250,
      dispersion_model = trend_model_observed,
      dataset = i
    )

    if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
    if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)
    if (!dir.exists(popeff_dir)) dir.create(popeff_dir, recursive = TRUE)

    if (!dir.exists(est_dir)) dir.create(est_dir, recursive = TRUE)
    if (!dir.exists(eval_dir)) dir.create(eval_dir, recursive = TRUE)

    # Resume per DATASET, not per directory. The previous guard asked whether
    # est_dir contained any files at all, which meant that once dataset 1 had
    # written its estimates the remaining nine were skipped -- a run from an
    # empty directory produced one dataset, not ten. Only overwrite = TRUE
    # forced the loop to complete, so the argument was doing the work the
    # guard was supposed to do. m5_evaluations is the last file written for a
    # dataset, so its presence means this dataset finished.
    done <- file.exists(
      file.path(eval_dir, paste0("m5_evaluations_", i, ".RDS"))
    )

    if (!done | overwrite) {
      # Save simulated data for later
      saveRDS(
        d$simdat,
        file = file.path(raw_dir, paste0("dataset_", i, ".RDS"))
      )
      saveRDS(
        d$combined_data,
        file = file.path(clean_dir, paste0("clean_dataset_", i, ".RDS"))
      )
      saveRDS(
        d$population_effects,
        file = file.path(popeff_dir, paste0("population_effects_", i, ".RDS"))
      )

      # Model 1 and 2 ##########################
      # Fitting naive and informed Neg-Binom model. The informed model has
      # a wider prior for the mean-dispersion fit as we use weighted estimates
      # of the loess regression.
      m1_m2_results <- m1_m2_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        weighted_loess = TRUE,
        CORES = cores
      )

      saveRDS(
        m1_m2_results$summaries_m1,
        file = file.path(est_dir, paste0("m1_estimates_", i, ".RDS"))
      )
      saveRDS(
        m1_m2_results$summaries_m2,
        file = file.path(est_dir, paste0("m2_estimates_", i, ".RDS"))
      )

      saveRDS(
        m1_m2_results$evaluations_m1,
        file = file.path(eval_dir, paste0("m1_evaluations_", i, ".RDS"))
      )
      saveRDS(
        m1_m2_results$evaluations_m2,
        file = file.path(eval_dir, paste0("m2_evaluations_", i, ".RDS"))
      )

      # Model 3 ################################################
      # A model for transformed counts.
      m3_results <- m3_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        CORES = cores
      )

      saveRDS(
        m3_results$summaries_m3,
        file = file.path(est_dir, paste0("m3_estimates_", i, ".RDS"))
      )
      saveRDS(
        m3_results$evaluations_m3,
        file = file.path(eval_dir, paste0("m3_evaluations_", i, ".RDS"))
      )

      # Model 4 and 5 ################################################
      # This is the Poisson model with observation-level random effects
      # Model 5 is the informed model (with priors).
      m4_m5_results <- m4_m5_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        CORES = cores
      )

      saveRDS(
        m4_m5_results$summaries_m4,
        file = file.path(est_dir, paste0("m4_estimates_", i, ".RDS"))
      )
      saveRDS(
        m4_m5_results$summaries_m5,
        file = file.path(est_dir, paste0("m5_estimates_", i, ".RDS"))
      )

      saveRDS(
        m4_m5_results$evaluations_m4,
        file = file.path(eval_dir, paste0("m4_evaluations_", i, ".RDS"))
      )
      saveRDS(
        m4_m5_results$evaluations_m5,
        file = file.path(eval_dir, paste0("m5_evaluations_", i, ".RDS"))
      )

      print(paste0("Simulation 1:", i, " is done."))
    }
  }
}


#' Simulation wrapper 2
#'
#' A wrapper for the simulation functions for scenario 2
#'
#' @param cores Number of cores to be used
#' @param seed Set the seed.
#' @param overwrite Should available results be overwritten?
sim_wrap2 <- function(cores, seed = 1, overwrite = FALSE) {
  set.seed(seed)

  .fw <- fast_workers_begin()
  on.exit(fast_workers_end(.fw), add = TRUE)

  raw_dir <- here::here("analysis/data/raw_data/simdata2/raw")
  clean_dir <- here::here("analysis/data/raw_data/simdata2/clean")
  popeff_dir <- here::here("analysis/data/raw_data/simdata2/popeffect")
  est_dir <- here::here("analysis/data/raw_data/estimates2")
  eval_dir <- here::here("analysis/data/raw_data/evaluations2")

  for (i in 1:10) {
    d <- simulate_datasets(
      nullgenes = 7500,
      condB_true = 1250,
      condB_time2_true = 1250,
      dispersion_model = trend_model_observed_low,
      dataset = i
    )

    if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
    if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)
    if (!dir.exists(popeff_dir)) dir.create(popeff_dir, recursive = TRUE)

    if (!dir.exists(est_dir)) dir.create(est_dir, recursive = TRUE)
    if (!dir.exists(eval_dir)) dir.create(eval_dir, recursive = TRUE)

    # Resume per DATASET, not per directory. The previous guard asked whether
    # est_dir contained any files at all, which meant that once dataset 1 had
    # written its estimates the remaining nine were skipped -- a run from an
    # empty directory produced one dataset, not ten. Only overwrite = TRUE
    # forced the loop to complete, so the argument was doing the work the
    # guard was supposed to do. m5_evaluations is the last file written for a
    # dataset, so its presence means this dataset finished.
    done <- file.exists(
      file.path(eval_dir, paste0("m5_evaluations_", i, ".RDS"))
    )

    if (!done | overwrite) {
      # Save simulated data for later
      saveRDS(
        d$simdat,
        file = file.path(raw_dir, paste0("dataset_", i, ".RDS"))
      )
      saveRDS(
        d$combined_data,
        file = file.path(clean_dir, paste0("clean_dataset_", i, ".RDS"))
      )
      saveRDS(
        d$population_effects,
        file = file.path(popeff_dir, paste0("population_effects_", i, ".RDS"))
      )

      # Model 1 and 2 ##########################
      # Fitting naive and informed Neg-Binom model. The informed model has
      # a wider prior for the mean-dispersion fit as we use weighted estimates
      # of the loess regression.
      m1_m2_results <- m1_m2_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        weighted_loess = TRUE,
        CORES = cores
      )

      saveRDS(
        m1_m2_results$summaries_m1,
        file = file.path(est_dir, paste0("m1_estimates_", i, ".RDS"))
      )
      saveRDS(
        m1_m2_results$summaries_m2,
        file = file.path(est_dir, paste0("m2_estimates_", i, ".RDS"))
      )

      saveRDS(
        m1_m2_results$evaluations_m1,
        file = file.path(eval_dir, paste0("m1_evaluations_", i, ".RDS"))
      )
      saveRDS(
        m1_m2_results$evaluations_m2,
        file = file.path(eval_dir, paste0("m2_evaluations_", i, ".RDS"))
      )

      # Model 3 ################################################
      # A model for transformed counts.
      m3_results <- m3_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        CORES = cores
      )

      saveRDS(
        m3_results$summaries_m3,
        file = file.path(est_dir, paste0("m3_estimates_", i, ".RDS"))
      )
      saveRDS(
        m3_results$evaluations_m3,
        file = file.path(eval_dir, paste0("m3_evaluations_", i, ".RDS"))
      )

      # Model 4 and 5 ################################################
      # This is the Poisson model with observation-level random effects
      # Model 5 is the informed model (with priors).
      m4_m5_results <- m4_m5_sim(
        d$combined_data,
        dataset = i,
        dofit = TRUE,
        CORES = cores
      )

      saveRDS(
        m4_m5_results$summaries_m4,
        file = file.path(est_dir, paste0("m4_estimates_", i, ".RDS"))
      )
      saveRDS(
        m4_m5_results$summaries_m5,
        file = file.path(est_dir, paste0("m5_estimates_", i, ".RDS"))
      )

      saveRDS(
        m4_m5_results$evaluations_m4,
        file = file.path(eval_dir, paste0("m4_evaluations_", i, ".RDS"))
      )
      saveRDS(
        m4_m5_results$evaluations_m5,
        file = file.path(eval_dir, paste0("m5_evaluations_", i, ".RDS"))
      )

      print(paste0("Simulation 1:", i, " is done."))
    }
  }
}
