
# Functions used in simulations and analysis #########################
#
#
#
# 01. filter_fun Filtering expression of simulated genes to remove
# "low expression" genes.
# 02. simulate_datasets A wrapper function to simplify simulation of datasets
# 03. sigma_summary A summary function for seqwrap - used in NB models
# 04. lmer_summary Summary function for seqwrap for lmer-based models (is singular)
# 05. poisson_summary Summary function for Poisson models
# 06. simulation functions - Wrapper functions for specific models used in
# loop over all data sets
# 07. extract_simulations Extract simulations from dataset-specific files.

# 01. filter_fun ##############################################################
# A function for filter by expression, some genes will have low expression
# due to sampling variability, these are removed in this function.
#
# counts, a count data frame. First column are target id.
# metadata, a metadat data frame containing group/time combinations
filter_fun <- function(counts, metadata) {

  ## Filter by expression
  # Combine all gene counts after filtering
  keep <- filterByExpr(
    counts[,-1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(metadata$condition, metadata$time))

  counts_filtered <-  counts[keep,]


  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(counts_filtered[,-1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(- group)

  # Combine all meta data
  metadata <- metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate( efflibsize = (lib.size * norm.factors) /
              median(lib.size * norm.factors),
            ln_efflibsize = log(efflibsize))

  return(list(counts = counts_filtered,
              metadata = metadata))


}


# 02. simulate_datasets ########################################################
# A (wrapper) function for simulating data sets
#
# nullgenes, the number of genes with null effects
# condB_true, number of genes with non-zero effects in baseline group diffs
# condB_timet2_true, number of genes with non-zero effects in interaction
simulate_datasets <- function(nullgenes = 7500,
                              condB_true = 1250,
                              condB_time2_true = 1250,
                              dispersion_model = NULL,
                              dataset) {



  ngenes <- nullgenes + condB_true + condB_time2_true

  ## Set fixed (population level) effects
  beta0 <- runif(nullgenes + condB_true + condB_time2_true,
                 min = 1.5,
                 max = 7)
  conditionB <-  c(rep(0, nullgenes),
                   runif(
                     condB_true,
                     min = 0.2,
                     max = 1) * sample(c(-1, 1),
                                       condB_true,
                                       prob = c(0.5, 0.5),
                                       replace = TRUE),
                   rep(0, condB_time2_true))

  timet2 <-  rnorm(nullgenes + condB_true + condB_time2_true,
                   0,
                   0.1)
  timet3 <-  rnorm(nullgenes + condB_true + condB_time2_true,
                   0,
                   0.2)

  conditionB_timet2 <- rep(0, nullgenes +
                             condB_true +
                             condB_time2_true)

  conditionB_timet3 <-  c(rep(0, nullgenes),
                          rep(0, condB_true),
                          runif(
                            condB_true,
                            min = 0.2,
                            max = 1) * sample(c(-1, 1),
                                              condB_true,
                                              prob = c(0.25, 0.75),
                                              replace = TRUE))

  # Simulate random effects
  # approximately based on observed data
  b0_values <- rlnorm(nullgenes + condB_true + condB_time2_true,
                      meanlog = -2.07,
                      sdlog = 1)



  # Set b1 and b2 distribution to ~small
  b1_values <- rep(0, nullgenes + condB_true + condB_time2_true)
  b2_values <- rep(0, nullgenes + condB_true + condB_time2_true)


  # Simulate data #
  simdat <- simcounts2(n1 = 32,
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
                       max_prop = 0.02)

  # Subdivide data sets into different sample sizes
  # small = 8 + 10
  # medium = 16 + 20 (similar to the observed)
  # large = 32 + 40

  metadata_small <- simdat$metadata |>
    filter(id %in% c(paste0("A", 1:8),
                     paste0("B", 1:10)))
  metadata_medium <- simdat$metadata |>
    filter(id %in% c(paste0("A", 1:16),
                     paste0("B", 1:20)))
  metadata_large <- simdat$metadata |>
    filter(id %in% c(paste0("A", 1:32),
                     paste0("B", 1:40)))


  counts_small <- simdat$counts |>
    select(gene, all_of(metadata_small$seq_sample_id))

  counts_medium <- simdat$counts |>
    select(gene, all_of(metadata_medium$seq_sample_id))


  counts_large <- simdat$counts |>
    select(gene, all_of(metadata_large$seq_sample_id))


  # Filter low expression genes, the resulting count tables contain
  # all genes that are to be used in simulations.

  combined_data <- list(small =  filter_fun(counts = counts_small,
                                            metadata = metadata_small),
                        medium =  filter_fun(counts = counts_medium,
                                             metadata = metadata_medium),
                        large = filter_fun(counts = counts_large,
                                           metadata = metadata_large))

  # Genes present in data sets after filtering
  genes_small <- combined_data[[1]]$counts$gene
  genes_medium <- combined_data[[2]]$counts$gene
  genes_large <- combined_data[[3]]$counts$gene

  filtered_genes <- data.frame(size = c(rep("small", length(genes_small)),
                                        rep("medium", length(genes_medium)),
                                        rep("large", length(genes_large))),
                               target = c(genes_small,
                                          genes_medium,
                                          genes_large)) |>
    expand_grid(term = c("conditionB", "timet3:conditionB"))


  # Save population effects
  population_effects <- bind_rows(

    data.frame(target = 1:ngenes,
               term = rep("conditionB", ngenes),
               population_effect = conditionB,
               dataset = dataset),
    data.frame(target = 1:ngenes,
               term = rep("timet3:conditionB", ngenes),
               population_effect = conditionB_timet3,
               dataset = dataset)) |>
    expand_grid(size = c("small", "medium", "large")) |>
    inner_join(filtered_genes)

  return(list(
    simdat = simdat,
    combined_data = combined_data,
    population_effects = population_effects
              ))



}



# 03. sigma_summary ###########################################################
# A summary function to return the dispersion parameter with SE on the log scale
# mean(predict(x, type = "link)) will give us the predicted log counts.
# We will put this in the eval fun to also get estimates of the parameters in
# the generic summary function.
sigma_summary <- function(x) {

  if(is.null(x$fit$convergence)) {
    conv <- 1
  } else {
    conv <- x$fit$convergence[[1]]
  }


  out <- data.frame(dispersion =  data.frame(summary(x$sdr))["betadisp",1],
                    dispersion.se = data.frame(summary(x$sdr))["betadisp",2],
                    log_mu = mean(predict(x, type = "link")),
                    convergence = conv,
                    pdHess = x$sdr$pdHess )
  return(out)

}


# 04. lmer_summary ############################################################
# A summary function for the lmer model of transformed counts
# it will return a the singularity diagnostics from lme4.
lmer_summary <- function(x) {

  out <- data.frame(isSingular = lme4::isSingular(x))

  return(out)

}

# 05. poisson_summary #########################################################
# A summary function for the Poisson models. This will use the convergence
# diagnostic in glmmTMB to indicate convergence. The pdHess indicator from the
# Hessian matrices
# see https://stackoverflow.com/questions/79110546/glmmtmb-convergence-messages
poisson_summary <- function(x) {

  if(is.null(x$fit$convergence)) {
    conv <- 1
  } else {
    conv <- x$fit$convergence[[1]]
  }

  out <- data.frame(convergence = conv,
                    pdHess = x$sdr$pdHess)
  return(out)


}


# 06. simulation functions ####################################################
# Model 1 and 2 simulation function
#
#
# Model 1 and 2 are the naive and informed negative binomial models
# weighted_loess, should a weighted loess be used for mean-dispersion estimates
# dataset, when used in a loop dataset indicate the index
# dofit, if false only data management is done

m1_m2_sim <- function(combined_data, dataset, dofit = TRUE,
                      weighted_loess = TRUE,
                      CORES = 2) {

  evaluations <- list()
  summaries <- list()
  evaluations2 <- list()
  summaries2 <- list()


  for(k1 in seq_along(combined_data)) {

    ms1  <- seqwrap_compose(
      data = combined_data[[k1]]$counts,                  # These are the filtered counts
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = sigma_summary,
      targetdata = NULL,
      arguments = list(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1|id),
        family = glmmTMB::nbinom2)
    )

    if (dofit) {
      ms1_results <- seqwrap(
        ms1,
        return_models = FALSE,
        # subset = 1:1100,
        cores = CORES)

      ms1_sum <- seqwrap_summarise(ms1_results)

      evaluations[[k1]] <- ms1_sum$evaluations |>
        mutate(model = "m1",
               datasets = dataset,
               size = names(combined_data)[k1])

      summaries[[k1]] <- ms1_sum$summaries  |>
        mutate(model = "m1",
               datasets = dataset,
               size = names(combined_data)[k1])



    }

    ## Model 2 ##
    # get successful targets
    targets <- evaluations[[k1]] |>
      filter(!is.na(dispersion.se)) |>
      distinct(target) |>
      pull(target)



    ## Fitting a model for the mean-dispersion relationship

    # Fit a trend to the dispersion data from m1
    # save the data in a convenient format. Using the log_mu_raw (average
    # observed counts) allows for a prior on dispersion also for genes with
    # unsuccessful fits in model 1. First we gather all dispersion data.

    dispersion_dat <- evaluations[[k1]] |>
      filter(target %in% targets)

    # Calculate the raw observed counts from the data
    raw_log_counts <- data.frame(
      target = as.character( combined_data[[k1]]$counts[, 1]),
      log_mu_raw = log(rowMeans( combined_data[[k1]]$counts[,-1]))
    )


    # Adding log raw counts to the dispersion df for modeling.
    #
    dispersion_dat <- dispersion_dat |>
      inner_join(raw_log_counts)

    ## Preliminary plot ##
    # dispersion_dat |>
    #  ggplot(aes(log_mu, log_mu_raw)) + geom_point()


    # Fit a loess model, using log_mu_raw as the predictor

    if(weighted_loess) {
      trend_model <- loess(dispersion ~ log_mu_raw,
                           data = dispersion_dat,
                           span = 0.7,
                           weights = 1/(dispersion.se^2))
    } else {
      trend_model <- loess(dispersion ~ log_mu_raw,
                           data = dispersion_dat,
                           span = 0.7)

    }



    # Display the model estimate mean-dispersion parameter
    ## Preliminary plot ##
    # dispersion_dat |>
    #   mutate(pred_dispersion = predict(trend_model)) |>
    #   ggplot(aes(log_mu_raw, dispersion)) +
    #   geom_point() +
    #   geom_line(aes(log_mu_raw, pred_dispersion), color = "red",
    #             linewidth = 2)
    #


    # Predict dispersion for each gene based on log raw counts
    # and combine into a prior.

    dispersion_prior <- data.frame(gene =  combined_data[[k1]]$counts[,1],
                                   pred = round(
                                     predict(trend_model,
                                             newdata = data.frame(
                                               log_mu_raw =  log(
                                                 rowMeans(
                                                   combined_data[[k1]]$counts[,-1]))
                                             )
                                     ), 3),
                                   s = round(
                                     trend_model$s, 3)
    ) |>
      mutate(prior = paste0("normal(", pred, ",", s, ")"))



    # Gather all distributions of estimates that can be used as prior information
    # in subsequent model. Fixed effects are centered on 0. We will use the SD
    # for priors


    ## Preliminary plot ##
    #  ms1_sum$summaries |>
    #  filter(target %in% targets) |>
    #  select(target, term, estimate) |>
    #  ggplot(aes(estimate)) + geom_density() +
    #  facet_wrap(~ term, scales = "free")

    estimate_distributions <- summaries[[k1]] |>
      filter(target %in% targets) |>
      select(target, term, estimate) |>

      summarise(.by = term,
                m = mean(estimate),
                s = sd(estimate)) |>
      filter(!term %in%  c("(Intercept)", "sd__(Intercept)"))

    # Extract the random effects distribution to fit a gamma distribution
    random_sd_estimate <- summaries[[k1]] |>
      filter(target %in% targets,
             term == "sd__(Intercept)") |>
      select(target, term, estimate) |>
      pull(estimate)

    # The gamma distribution is parameterized using a shape and a rate
    # parameter. It looks like this prior will lead to a push towards 0,
    # consider adding a constant to push away from zero...
    # TODO this may needs testing.
    mean_sd <- mean(random_sd_estimate)
    var_sd <- var(random_sd_estimate)
    shape_param <- 2 # mean_sd^2 / var_sd





    # Here we prepare priors for the fixed effects, in this version all fixed
    # effects, except the intercept, will have regularizing priors corresponding to
    # the distributions of effects seen in the naive models.
    Priors_df <- bind_rows(
      data.frame(prior = paste0("normal(0,",round(estimate_distributions$s,2 ), ")"),
                 class = rep("fixef", 5),
                 coef = estimate_distributions$term),
      data.frame(prior = paste0(
        "gamma(",
        round(mean_sd ,2),
        ",",
        2,
        ")"),
        class = "ranef",
        coef = "id")
    )

    # We want to use the mean-dispersion relationship to add a prior for the
    # dispersion parameter. This means that we need a gene specific prior

    Priors_list <- list()
    for( j in 1:nrow( combined_data[[k1]]$counts)) {


      df <- bind_rows(Priors_df,
                      data.frame(
                        prior =
                          dispersion_prior[dispersion_prior$gene ==
                                             combined_data[[k1]]$counts[j,1],4],
                        class = "fixef_disp",
                        coef = "1"
                      )
      )

      Priors_list[[j]] <- df



    }


    # seqwrap accepts target-wise data frames as a list,
    # this makes it easier to specify target-specific priors.

    # Here we specify priors based the results from ms1

    ms2 <- seqwrap_compose(
      data =  combined_data[[k1]]$counts,
      metadata =  combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = sigma_summary,
      targetdata =  Priors_list,
      arguments = alist(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1|id),
        family = glmmTMB::nbinom2,
        priors = data.frame(
          prior = prior,
          class = class,
          coef = coef) ))

    ms2_results <- seqwrap(
      ms2,
      return_models = FALSE,
      #  subset = 1:50,
      cores = CORES)


    ms2_sum <- seqwrap_summarise(ms2_results)


    evaluations2[[k1]] <- ms2_sum$evaluations |>
      mutate(model = "m2",
             datasets = dataset,
             size = names(combined_data)[k1])

    summaries2[[k1]] <- ms2_sum$summaries  |>
      mutate(model = "m2",
             datasets = dataset,
             size = names(combined_data)[k1])




  }

 return(list(
   evaluations_m1 = bind_rows(evaluations),
   evaluations_m2 = bind_rows(evaluations2),
   summaries_m1 = bind_rows(summaries),
   summaries_m2 = bind_rows(summaries2)
 ))


}


# This is the model for transformed counts data
# the function performs transformation and fits models over a combined_data
# object.
m3_sim <- function(combined_data, dataset, dofit = TRUE, CORES = 2) {

  evaluations <- list()
  summaries <- list()

  for(k3 in seq_along(combined_data)) {
    # The VST transformation

    # Safe check of any NA's in the data
    safe_counts <- combined_data[[k3]]$counts[
      complete.cases(combined_data[[k3]]$counts[,-1]),-1]


    dds <- DESeqDataSetFromMatrix(countData = safe_counts,
                                  colData = combined_data[[k3]]$metadata,
                                  design = ~ time * condition)

    dds <- DESeq(dds, quiet = TRUE)
    vst_mat <- assay(varianceStabilizingTransformation(dds,
                                                       blind = FALSE,
                                                       fitType = "parametric"))
    vst_dat <- cbind(data.frame(gene =  combined_data[[k3]]$counts[
      complete.cases(combined_data[[k3]]$counts[,-1]),1],
      as.data.frame(vst_mat)))



    ms3 <-  seqwrap_compose(
      data = vst_dat,
      metadata =  combined_data[[k3]]$metadata,
      samplename = "seq_sample_id",
      modelfun = lmerTest::lmer,
      eval_fun = lmer_summary,
      arguments = list(
        formula = y ~ time * condition + (1|id))
    )

    if (dofit) {
      ms3_results <- seqwrap(
        ms3,
        return_models = FALSE,
        #   subset = 1:50,
        cores = CORES)



      ms3_sum <- seqwrap_summarise(ms3_results)

      evaluations[[k3]] <- ms3_sum$evaluations |>
        mutate(model = "m3",
               datasets = dataset,
               size = names(combined_data)[k3])
      summaries[[k3]] <- ms3_sum$summaries|>
        mutate(model = "m3",
               datasets = dataset,
               size = names(combined_data)[k3])



    }


  }

  evaluations <- bind_rows(evaluations)
  summaries <- bind_rows(summaries)

  return(list(
    evaluations_m3 = bind_rows(evaluations),
    summaries_m3 = bind_rows(summaries)
  ))

}

# This function does the Poisson models with observation level random effects
# both in a naive and informed version.
m4_m5_sim <- function(combined_data, dataset, dofit = TRUE,
                      CORES = 2) {

  evaluations <- list()
  summaries <- list()
  evaluations2 <- list()
  summaries2 <- list()


  for(k1 in seq_along(combined_data)) {


    ms1  <- seqwrap_compose(
      data = combined_data[[k1]]$counts,                  # These are the filtered counts
      metadata = combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = poisson_summary,
      targetdata = NULL,
      arguments = list(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1|id) + (1|seq_sample_id),
        family = stats::poisson)
    )

    if (dofit) {
      ms1_results <- seqwrap(
        ms1,
        return_models = FALSE,
        # subset = 1:1100,
        cores = CORES)

      ms1_sum <- seqwrap_summarise(ms1_results)

      evaluations[[k1]] <- ms1_sum$evaluations |>
        mutate(model = "m4",
               datasets = dataset,
               size = names(combined_data)[k1])

      summaries[[k1]] <- ms1_sum$summaries  |>
        mutate(model = "m4",
               datasets = dataset,
               size = names(combined_data)[k1])


    }

    ## Model 2 ##
    # get successful targets
    targets <- evaluations[[k1]] |>
      filter(convergence == 0) |>
      distinct(target) |>
      pull(target)




    ## Preliminary plot ##
    #  ms1_sum$summaries |>
    #  filter(target %in% targets) |>
    #  select(target, term, estimate) |>
    #  ggplot(aes(estimate)) + geom_density() +
    #  facet_wrap(~ term, scales = "free")

    estimate_distributions <- summaries[[k1]] |>
      filter(target %in% targets) |>
      select(target, term, group, estimate) |>

      summarise(.by = c(term,group),
                m = mean(estimate),
                s = sd(estimate)) |>
      filter(!term %in%  c("(Intercept)", "sd__(Intercept)"))

    # Extract the random effects distribution to fit a gamma distribution
    # on the participant level intercepts
    random_sd_estimate <- summaries[[k1]] |>
      filter(target %in% targets,
             term == "sd__(Intercept)",
             group == "id") |>
      select(target, term, estimate) |>
      pull(estimate)

    random_sd_estimate_obs <- summaries[[k1]] |>
      filter(target %in% targets,
             term == "sd__(Intercept)",
             group == "seq_sample_id") |>
      select(target, term, estimate) |>
      pull(estimate)



    # The gamma distribution is parameterized using a shape and a rate
    # parameter. It looks like this prior will lead to a push towards 0,
    # consider adding a constant to push away from zero...
    # TODO this may needs testing.
    mean_sd <- mean(random_sd_estimate)
    var_sd <- var(random_sd_estimate)
    shape_param <- 2 # mean_sd^2 / var_sd
    # Observation level
    mean_sd_obs <- mean(random_sd_estimate_obs)
    var_sd_obs <- var(random_sd_estimate_obs)
    shape_param_obs <- 2 # mean_sd^2 / var_sd




    # Here we prepare priors for the fixed effects, in this version all fixed
    # effects, except the intercept, will have regularizing priors corresponding to
    # the distributions of effects seen in the naive models.

    # The priors for the random effects are added, glmmTMB accepts coefficient
    # index.
    Priors_df <- bind_rows(
      data.frame(prior = paste0("normal(0,",round(estimate_distributions$s,2 ), ")"),
                 class = rep("fixef", 5),
                 coef = estimate_distributions$term),
      data.frame(prior = paste0(
        "gamma(",
        c(round(mean_sd ,2),round(mean_sd_obs ,2)),
        ",",
        2,
        ")"),
        class = "ranef",
        coef = c("1", "2"))
    )

    # We want to use the mean-dispersion relationship to add a prior for the
    # dispersion parameter. This means that we need a gene specific prior

    Priors_list <- list()
    for( j in 1:nrow( combined_data[[k1]]$counts)) {

      Priors_list[[j]] <- Priors_df

    }


    # seqwrap accepts target-wise data frames as a list,
    # this makes it easier to specify target-specific priors.

    # Here we specify priors based the results from ms1

    ms2 <- seqwrap_compose(
      data =  combined_data[[k1]]$counts,
      metadata =  combined_data[[k1]]$metadata,
      samplename = "seq_sample_id",
      modelfun = glmmTMB::glmmTMB,
      eval_fun = poisson_summary,
      targetdata =  Priors_list,
      arguments = alist(
        formula = y ~ time * condition + offset(ln_efflibsize) + (1|id) + (1|seq_sample_id),
        family = stats::poisson,
        priors = data.frame(
          prior = prior,
          class = class,
          coef = coef) ))

    ms2_results <- seqwrap(
      ms2,
      return_models = FALSE,
      #  subset = 1:50,
      cores = CORES)


    ms2_sum <- seqwrap_summarise(ms2_results)


    evaluations2[[k1]] <- ms2_sum$evaluations |>
      mutate(model = "m5",
             datasets = dataset,
             size = names(combined_data)[k1])

    summaries2[[k1]] <- ms2_sum$summaries  |>
      mutate(model = "m5",
             datasets = dataset,
             size = names(combined_data)[k1])




  }

  return(list(
    evaluations_m4 = bind_rows(evaluations),
    evaluations_m5 = bind_rows(evaluations2),
    summaries_m4 = bind_rows(summaries),
    summaries_m5 = bind_rows(summaries2)
  ))


}


# 07. extract_simulations #####################################################

extract_simulations <- function(evaluations_path = "data_sim/evaluations",
                                estimates_path = "data_sim/estimates",
                                populationeffects_path = "data_sim/simdata/popeffect",
                                disp_scenario = "s1") {

  # Evaluations

  # The evaluations can be combined despite having different shape. Non-
  # available columns will be NA.

  eval_files <- list.files(evaluations_path)
  evaluations <- list()
  for(i in seq_along(eval_files)) {

    evaluations[[i]] <- readRDS(paste0(evaluations_path, "/",
                                       eval_files[i])) |>
      mutate(file = eval_files[i],
             disp_scenario = disp_scenario)

  }

  evaluations <- bind_rows(evaluations)




  # Estimates
  est_files <- list.files(estimates_path)
  estimates <- list()
  for(i in seq_along(est_files)) {

    estimates[[i]] <- readRDS(paste0(estimates_path, "/",
                                     est_files[i])) |>
      mutate(file = est_files[i],
             disp_scenario = disp_scenario)

  }

  estimates <- bind_rows(estimates)

  # Population effects
  pop_files <- list.files(populationeffects_path)
  popeffects <- list()
  for(i in seq_along(pop_files)) {

    popeffects[[i]] <- readRDS(paste0(populationeffects_path, "/",
                                      pop_files[i])) |>
      mutate(file = pop_files[i],
             disp_scenario = disp_scenario)

  }

  popeffects <- bind_rows(popeffects)

  return(list(populationeffects = popeffects,
              estimates = estimates,
              evaluations = evaluations))

}





