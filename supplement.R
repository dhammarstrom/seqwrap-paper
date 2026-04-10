#| label: check-packages-dependencies
#| echo: false
#| message: false
#| warning: false
#| output: false

# Being able to re-produce this document requires a set of packages. 
source("R/check-packages.R")

# Output needed for the figures are created in data.prep.R and 
# m1-m5-pillon.R
source("R/data-prep.R")
## Only needed if output data is missong
if(!file.exists("data-out/pillon-models.RDS")) {
  source("R/m1-m5-pillon-data.R")
}


#| eval: false
# 
# # Installing from GitHub requires remotes
# # install.packages("remotes")
# 
# remotes::install_github("trainome/seqwrap")
# 

#| message: false
#| warning: false
#| labels: packages-and-data

# Load packages
library(tidyverse)
library(seqwrap)
library(edgeR)
library(gt)

# Load data 
d <- dungan_counts

metadata <- dungan_counts$metadata
counts <- dungan_counts$countdata


#| label: tbl-meta-dungan
#| echo: false
#| message: false
#| warning: false
#| tbl-cap: Selected rows from the Dungan et al. metadata data frame.


metadata |> 
  slice_head(n = 2, by = c(treatment, surgery)) |> 
  gt() |> 
  opt_table_font(size = 12)

# Save seq_sample_id for display
samples <- metadata |> 
  slice_head(n = 2, by = c(treatment, surgery)) |> 
  pull(seq_sample_id)


#| label: tbl-dungan-counts
#| echo: false
#| message: false
#| warning: false
#| tbl-cap: Selected rows from the Dungan et al. counts data frame.

counts |> 
  slice_head(n = 10) |> 
  dplyr::select(gene_name, all_of(samples)) |> 
  gt() |> 
  opt_table_font(size = 12)




# Find targets with more than a unique row in the data
rm_targets <- counts |> 
  summarise(.by = gene_name, 
            n = n()) |> 
  filter(n > 1) |> 
  pull(gene_name)

# Filter out non-uniqe rows
counts <- counts |> 
  filter(!(gene_name %in% rm_targets))


#| echo: false

# Count the number of targets
n_genes <- nrow(counts)




# Use EdgeR to calculate the TMM
y <- edgeR::DGEList(counts[,-1])
y <- edgeR::calcNormFactors(y)

# Combine the data into metadata
metadata <- metadata |> 
  inner_join(
    y$samples |> 
  tibble::rownames_to_column("seq_sample_id")
  ) |> 
   mutate(
     efflibsize = (lib.size * norm.factors) / median(lib.size * norm.factors),
     efflibsize = log(efflibsize)
   ) 





m1.1 <- seqwrap_compose(
  modelfun = glmmTMB::glmmTMB,
  data = counts, 
  metadata = metadata, 
  arguments = list(formula = y ~ treatment * surgery + offset(efflibsize), 
                family = glmmTMB::nbinom2),
  samplename = "seq_sample_id"
  )




m1.1_temp <- seqwrap(
  m1.1, 
  cores = 1, 
  subset = 1:24, 
  return_models = TRUE, 
  verbose = FALSE
)




m1.1_temp_sum <- seqwrap_summarise(m1.1_temp, verbose = FALSE)


#| label: tbl-dungan-sum-temp
#| echo: false
#| message: false
#| warning: false
#| tbl-cap: The first ten rows from summarised preliminary models (model summaries).


m1.1_temp_sum$summaries |> 
  slice_head(n = 10) |> 
  gt() |> 
  fmt_auto() |> 
  opt_table_font(size = 10)


#| label: tbl-dungan-eval-temp
#| echo: false
#| message: false
#| warning: false
#| tbl-cap: The first ten rows from summarised preliminary models (model evaluations).
#|   Values are p-values from the DHARMa package (see [here for details](https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html)).


m1.1_temp_sum$evaluations |> 
  slice_head(n = 10) |> 
  gt() |> 
  opt_table_font(size = 12)




dispersion_fun <- function(x) {
  
    out <- data.frame(
      dispersion =  data.frame(summary(x$sdr))["betadisp",1],
      dispersion.se = data.frame(summary(x$sdr))["betadisp",2],
      log_mu = log(mean(x$frame$y))
      )
  
    return(out)
  
}


dispersion_fun(m1.1_temp@models[[1]])


log(sigma(m1.1_temp@models[[1]]))




m1.2_results <- seqwrap(
  m1.1,
  eval_fun = dispersion_fun,
  cores = 12, 
  return_models = FALSE,
  verbose = FALSE
)





warnings <- m1.2_results@errors|> 
  # filter the non-null elements
  filter(map_lgl(warnings_fit, ~ !is.null(.x[[1]]))) |> 
  # save in object
  pull(warnings_fit)

warnings[1]




m1.2_sum <- seqwrap_summarise(m1.2_results, verbose = FALSE)


#| label: fig-dungan-estimates
#| echo: false
#| message: false
#| warning: false
#| fig-cap: Estimates from modelling the Dungan data. The dispersion trend is visualized
#|   using the default `ggplot`2::geom_smooth`.

library(cowplot)
library(ggtext)


p1 <- m1.2_sum$evaluations |> 
  ggplot(aes(log_mu, dispersion)) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE, color = "orange") + 
  theme_classic() +
  labs(x = "Average raw counts log(&mu;) ", 
       y = "Dispersion log(&theta;)") + 
  theme(axis.title.x = element_markdown(), 
        axis.title.y = element_markdown())
  

p2 <- m1.2_sum$summaries |> 
  
  ggplot(aes(estimate)) + 
  geom_density() + 
  facet_wrap(~ term, scales = "free") +
  theme_classic() + 
  labs(x = "Estimate", y = "")


plot_grid(
  plot_grid(NULL, p1, NULL, rel_heights = c(0.2, 1, 0.2), ncol = 1),
            p2, 
          ncol = 2, 
          rel_widths = c(1, 1.5))




# Summarise parameters from the main output
param_sum <- m1.2_sum$summaries |> 
  summarise(.by = term, 
            m = mean(estimate), 
            s = sd(estimate)) 

# Fit a model to the log dispersion values
trend_model <- loess(dispersion ~ log_mu,
                         data = m1.2_sum$evaluations,
                         span = 0.7)


# Predict dispersion for each gene based on log raw counts
# and combine into a prior.
dispersion_prior <- data.frame(gene =  counts[,1],
                               pred = round(
                                 predict(trend_model,
                                         newdata = data.frame(
                                           log_mu =  log(rowMeans(counts[,-1]))
                                           )
                                   ), 3),
                                 s = round(
                                   trend_model$s, 3)
  ) |>
    mutate(prior = paste0("normal(", pred, ",", s, ")"))

# Combine parameters in a data frame...
Priors_df <- data.frame(prior = paste0("normal(",
                                       round(param_sum$m, 2),
                                       ",",
                                       round(param_sum$s,2 ), 
                                       ")"),
                        class = rep("fixef", 4),
                        coef = param_sum$term)


# ... extract a prediction for the mean dispersion parameter per target
# and combine with other parameters  
Priors_list <- list()
for( j in 1:nrow( counts )) {


    df <- bind_rows(Priors_df,
                    data.frame(
                      prior =
                        dispersion_prior[dispersion_prior$gene ==
                                           counts[j,1],4],
                      class = "fixef_disp",
                      coef = "1"
                    )
    )

    Priors_list[[j]] <- df

}





m1.3 <- seqwrap_compose(
  modelfun = glmmTMB::glmmTMB,
  data = counts, 
  metadata = metadata, 
  eval_fun = dispersion_fun, 
  arguments =  alist(
    formula = y ~ treatment * surgery + offset(efflibsize), 
    family = glmmTMB::nbinom2,
      priors = data.frame(
        prior = prior,
        class = class,
        coef = coef)), 
  targetdata = Priors_list
)


m1.3_temp <- seqwrap(
  m1.3, 
  subset = 1:24, 
  return_models = TRUE, 
  verbose = FALSE
)


summary(m1.3_temp@models[[1]])




m1.3_results <- seqwrap(
  m1.3, 
  cores = 12,
  return_models = FALSE, 
  verbose = FALSE
)

m1.3_sum <- seqwrap_summarise(m1.3_results, verbose = FALSE)


#| label: fig-dungan-final-estimates
#| echo: false
#| message: false
#| warning: false
#| fig-cap: A visualization of results from the regularized model of the Dungan et al.
#|   data.



m1.3_sum$summaries |> 

  filter(term == "treatmentSenolytic:surgeryOverload") |> 
  mutate(fdr = p.adjust(p.value, method = "fdr"), 
         sig = if_else(fdr < 0.05, "significant", "non-significant"), 
         SYMBOL = if_else(sig == "significant", target, "")) |> 
  
  ggplot(aes(estimate, -log2(p.value), 
             color = sig)) + 
  geom_point(alpha = 0.4) + 
  
  geom_text(aes(label = SYMBOL), 
      show.legend = FALSE,
            position = position_nudge(y = 1)) +
  
  theme_classic() + 
  labs(x = "Estimate (log-scaled)", 
       y = "-log<sub>2</sub>(p-value)", 
       color = "FDR < 0.05", 
       subtitle = "Overload:Senolytic treatment") + 
  theme(axis.title.y = element_markdown()) 

  
  



#| label: comp-data
#| echo: true
#| message: false
#| warning: false

# Need installation from github: https://github.com/stop-pre16/lmerSeq
library(lmerSeq) 
library(DESeq2)
library(edgeR)

library(dplyr)
library(tidyverse)
library(seqwrap)
library(purrr)
library(glmmSeq)
library(ggplot2)

# Load the data from seqwrap
dat <- seqwrap::pillon_counts


# Filter low expressed genes and combine after filtering
keep <- filterByExpr(
  dat$countdata[,-1],
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7,
  group = paste(dat$metadata$group, dat$metadata$time))

countdat <-  dat$countdata[keep,]

metadat <- dat$metadata


# Check if sequence_ID matches.
all(colnames(countdat)[-1] == metadat$seq_sample_id)


# Transform the count data using DESeq2's VST
dds <- DESeqDataSetFromMatrix(countData = countdat[,-1],
                              colData = metadat,
                              design = ~ time * group)

dds <- DESeq(dds)
vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
vst_expr <- assay(vsd.fixed)





# Add as gene_ids to the VST-transformed counts
rownames(vst_expr) <- countdat$geneid

# subset to only the first fifteen genes
vst_expr <- vst_expr[1:15,]


#| label: lmerseq
#| echo: true
#| message: false
#| warning: false

# fit the model using lmerseq
fit.lmerSeq <- lmerSeq.fit(form = ~ time * group  + (1|id),
                           expr_mat = vst_expr,
                           sample_data = metadat,
                           parallel = FALSE,
                           REML = TRUE)



# Extract coefficients from the models 
coefs <- c("timepost", "timerec", "groupT2D", "timepost:groupT2D", "timerec:groupT2D")

# Extract all coefficients from lmerseq models
lmerseq_summary <- coefs |> 
  map_df(~ lmerSeq.summary(
    lmerSeq_results = fit.lmerSeq,
    coefficient = .x,
    p_adj_method = "BH",
    ddf = "Satterthwaite",
    sort_results = FALSE
  )$summary_table |> 
    mutate(Coefficient = .x)
  )




#| label: seqwrap
#| echo: true
#| message: false
#| warnings: false

# make data suitable for seqwrap as it accepts only data frames or lists
vst_expr_df <- data.frame(vst_expr) %>%
  tibble::rownames_to_column(var = "target")

# Creating a evaluation function for finding singular fits
singular_eval <- function(x) {
  out <- data.frame(isSingular = lme4::isSingular(x))
  return(out)
}



# first validate input using seqwrap compose
seqwrap_val <- seqwrap_compose(
  data = vst_expr_df,
  metadata = metadat,
  samplename = "seq_sample_id",
  modelfun = lmerTest::lmer,
  arguments = list(
    formula = y ~ time * group  + (1|id)), 
  eval_fun = singular_eval
)


seqwrap_model <- seqwrap(seqwrap_val,
                        return_models = FALSE,
                        cores = 1, 
                        verbose = FALSE)



# Check for singular fits and list only targets without 
targets <- seqwrap_summarise(seqwrap_model, 
                             verbose = FALSE)$evaluations |> 
  filter(!isSingular) |> 
  pull(target)


seqwrap_summary_df <- seqwrap_summarise(seqwrap_model, 
                                        verbose = FALSE)$summaries |> 
  filter(target %in% targets)



#| label: fig-comp-visualisation
#| echo: false
#| fig-cap: Comparing estimates from `lmerseq` and `seqwrap`
#| fig-height: 8


# Select
seqwrap_df <- seqwrap_summary_df %>%
  subset(!term %in% c("(Intercept)", "sd__(Intercept)", "sd__Observation")) %>%
  dplyr::select("target",  "estimate", "p.value", "term" ) %>%
  mutate(source = "seqwrap")


lmerseq_df <- lmerseq_summary |>
  dplyr::select(gene, Estimate, p_val_raw, Coefficient) |>
  mutate(source = "lmerseq") |>
  #rename to match the columns from seqwrap summary
  dplyr::rename(target = gene,
         estimate = Estimate,
         p.value = p_val_raw,
         term = Coefficient)



p1 <- bind_rows(seqwrap_df, 
          lmerseq_df) |> 
  pivot_wider(names_from = source, 
              values_from = c(estimate, p.value)) |> 
 ggplot(aes(x = estimate_seqwrap, y = estimate_lmerseq, color = target)) +
   geom_point(size = 2) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
   facet_wrap(~ term) +
   labs(subtitle = "Comparison of Estimates", x = "seqwrap Estimate", y = "lmerseq Estimate", 
        color = "") +
   theme_classic()



p2 <-  bind_rows(seqwrap_df, 
          lmerseq_df) |> 
  pivot_wider(names_from = source, 
              values_from = c(estimate, p.value)) |> 
   ggplot(aes(x = p.value_seqwrap, y = p.value_lmerseq, color = target)) +
     geom_point(size = 2) +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
     facet_wrap(~ term) +
   labs(subtitle = "Comparison of p-values", x = "seqwrap p-values", y = "lmerseq p-values") +
     theme_classic() + 
  theme(legend.position = "none")


plot_grid(p1, p2, ncol = 1)


#| eval: false
# 
# ### Code chunk not evaluated ###
# 
# 
# # Load packages
# library(seqwrap)
# library(tidyverse)
# library(edgeR)
# 
# # Loading the data. The Pillon 2022 data set is part of the seqwrap package as `pillon_counts`.
# dat <- seqwrap::pillon_counts
# 
# # Filter low expressed genes and combine after filtering
# keep <- filterByExpr(
#   dat$countdata[,-1],
#   min.count = 10,
#   min.total.count = 15,
#   large.n = 10,
#   min.prop = 0.7,
#   group = paste(dat$metadata$group, dat$metadata$time))
# 
# 
# 
# countdat <-  dat$countdata[keep,]
# metadat <- dat$metadata
# 
# 
# # Use EdgeR to calculate the TMM
# y <- edgeR::DGEList(countdat[,-1])
# y <- edgeR::calcNormFactors(y)
# 
# # Store effective library sizes 
# libsize <- y$samples |>
#     rownames_to_column(var = "seq_sample_id") |>
#     select(- group)
# 
# # Combine all meta data
# metadat <- dat$metadata |>
#   inner_join(libsize, by = "seq_sample_id") |>
#   mutate(group = factor(group, levels = c("NGT", "T2D")),
#          time = factor(time, levels = c("basal", "post", "rec")),
#          efflibsize = (lib.size * norm.factors)/median(lib.size * norm.factors),
#          ln_efflibsize = log(efflibsize))
# 
# # Setting up the evaluation function
# sigma_summary2 <- function(x) {
# 
#   if(is.null(x$fit$convergence)) {
#     conv <- 1
#   } else {
#     conv <- x$fit$convergence[[1]]
#   }
# 
# 
#   ### Simulating scaled resid ###
#   # This creates simulated residuals for test of uniformity and dispersion.
#   # See
#   # https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
#   # for details.
#   residObj <- DHARMa::simulateResiduals(x, plot = FALSE)
# 
#   # Saving tests from DHARMa
#   unif <- DHARMa::testUniformity(residObj, plot = FALSE)
#   disp <- DHARMa::testDispersion(residObj, plot = FALSE)
# 
# 
#   # Saving SD of Obs Level Rand Eff in the Poisson model
# 
#   tidy_coef <- broom.mixed::tidy(x)
# 
#   if (any(tidy_coef$group %in% "seq_sample_id")) {
#     olre.sd <- tidy_coef |>
#       dplyr::filter(group == "seq_sample_id") |>
#       dplyr::pull(estimate)
#   } else { olre.sd <- NA }
# 
# 
# 
#   # Combine all in a data frame
#   out <- data.frame(dispersion =  data.frame(summary(x$sdr))["betadisp",1],
#                     dispersion.se = data.frame(summary(x$sdr))["betadisp",2],
#                     olre.sd = olre.sd,
#                     log_mu = mean(predict(x, type = "link")),
#                     convergence = conv,
#                     pdHess = x$sdr$pdHess,
#                     aic = AIC(x),
#                     unif.p = unif$p.value,
#                     unif.stat = unif$statistic,
#                     disp.p = disp$p.value,
#                     disp.stat = disp$statistic)
#   return(out)
# 
# }
# 
# 
# 
# 
# # Setting up the model
# m1 <- seqwrap_compose(
#   data = countdat,
#   metadata = metadat,
#   samplename = "seq_sample_id",
#   modelfun = glmmTMB::glmmTMB,
#   eval_fun = sigma_summary2,
#   arguments = list(
#     formula = y ~ time * group + offset(ln_efflibsize) + (1|id),
#     family = glmmTMB::nbinom2)
#   )
# 
# m1_results <- seqwrap(m1,
#                       return_models = FALSE,
#                       cores = parallel::detectCores())
# 

#| eval: false
# 
# ### Code chunk not evaluated ###
# 
# m2 <- seqwrap_compose(
#   data = countdat,
#   metadata = metadat,
#   samplename = "seq_sample_id",
#   modelfun = glmmTMB::glmmTMB,
#   eval_fun = sigma_summary2,
#   targetdata = NULL,
#   arguments = list(
#     formula = y ~ time * group + 
#       offset(ln_efflibsize) + 
#       (1|id) +
#       (1|seq_sample_id),
#       family = stats::poisson)
#   )
# 
# 
# m2_results <- seqwrap(
#       m2,
#       return_models = FALSE,
#       verbose = FALSE,
#       subset = 1:16,
#       cores = parallel::detectCores())
# 

#| eval: false
# 
# ### Code not evaluated ###
# 
# ## Summarise results from the non-informed model
# m1_sums  <- seqwrap_summarise(m1_results)
# 
# 
# # Prepare priors for the regularized version
# 
#   ## Fitting a model for the mean-dispersion relationship
# 
#   # Fit a trend to the dispersion data from m1
#   # save the data in a convenient format. Using the log_mu_raw (average
#   # observed counts) allows for a prior on dispersion also for genes with
#   # unsuccessful fits in model 1. First we gather all dispersion data.
# 
#   dispersion_dat <- m1_sums$evaluation |>
#     filter(target %in% targets)
# 
#   # Calculate the raw observed counts from the data
#   raw_log_counts <- data.frame(
#     target = as.character( countdat[, 1]),
#     log_mu_raw = log(rowMeans(countdat[,-1]))
#   )
# 
# 
#   # Adding log raw counts to the dispersion df for modeling.
#   dispersion_dat <- dispersion_dat |>
#     inner_join(raw_log_counts)
# 
# 
# 
#   # The SE of the dispersion estimate is used as weights in the trended model
#     trend_model <- loess(dispersion ~ log_mu_raw,
#                          data = dispersion_dat,
#                          span = 0.7,
#                          weights = 1/(dispersion.se^2))
# 
# 
#   # Predict dispersion for each gene based on log raw counts
#   # and combine into a prior.
# 
#   dispersion_prior <- data.frame(gene =  countdat[,1],
#                                  pred = round(
#                                    predict(trend_model,
#                                            newdata = data.frame(
#                                              log_mu_raw =  log(
#                                                rowMeans(
#                                                  countdat[,-1]))
#                                            )
#                                    ), 3),
#                                  s = round(
#                                    trend_model$s, 3)
#   ) |>
#     mutate(prior = paste0("normal(", pred, ",", s, ")"))
# 
# 
# 
#   # Gather all distributions of estimates that can be used as prior information
#   # in subsequent model.
# 
# 
#   ## Preliminary plot ##
#   #  ms1_sum$summaries |>
#   #  filter(target %in% targets) |>
#   #  select(target, term, estimate) |>
#   #  ggplot(aes(estimate)) + geom_density() +
#   #  facet_wrap(~ term, scales = "free")
# 
#   estimate_distributions <- m1_sums$summaries |>
#     filter(target %in% targets) |>
#     select(target, term, estimate) |>
# 
#     summarise(.by = term,
#               m = mean(estimate),
#               s = sd(estimate)) |>
#     filter(!term %in%  c("(Intercept)", "sd__(Intercept)"))
# 
#   # Extract the random effects distribution to fit a gamma distribution
#   random_sd_estimate <- m1_sums$summaries |>
#     filter(target %in% targets,
#            term == "sd__(Intercept)") |>
#     select(target, term, estimate) |>
#     pull(estimate)
# 
#   # The gamma distribution is parameterized using a shape and a rate
#   # parameter. It looks like this prior will lead to a push towards 0,
#   # we are adding a constant to push away from zero...
#   mean_sd <- mean(random_sd_estimate)
#   var_sd <- var(random_sd_estimate)
#   shape_param <- 2 # mean_sd^2 / var_sd
# 
# 
# 
# 
# 
#   # Here we prepare priors for the fixed effects, in this version all fixed
#   # effects, except the intercept, will have regularizing priors 
#   # corresponding to the distributions of effects seen in the non-informed
#   # models.
#   Priors_df <- bind_rows(
#     data.frame(prior = paste0("normal(",
#                               round(estimate_distributions$m,2 ), 
#                               ", ",
#                               ,round(estimate_distributions$s,2 ), ")"),
#                class = rep("fixef", 5),
#                coef = estimate_distributions$term),
#     data.frame(prior = paste0(
#       "gamma(",
#       round(mean_sd ,2),
#       ",",
#       2,
#       ")"),
#       class = "ranef",
#       coef = "id")
#   )
# 
#   # We want to use the mean-dispersion relationship to add a prior for the
#   # dispersion parameter. This means that we need a target specific prior
# 
#   Priors_list <- list()
#   for( j in 1:nrow( countdat)) {
# 
# 
#     df <- bind_rows(Priors_df,
#                     data.frame(
#                       prior =
#                         dispersion_prior[dispersion_prior$gene ==
#                                            countdat[j,1],4],
#                       class = "fixef_disp",
#                       coef = "1"
#                     )
#     )
# 
#     Priors_list[[j]] <- df
# 
# 
# 
#   }
# 
# # Fit regularized model 
#   
#   m1.reg <- seqwrap_compose(
#     data =  countdat,
#     metadata =  metadat,
#     samplename = "seq_sample_id",
#     modelfun = glmmTMB::glmmTMB,
#     eval_fun = sigma_summary2,
#     targetdata =  Priors_list,
#     arguments = alist(
#       formula =  y ~ time * group + offset(ln_efflibsize) + (1|id),
#       family = glmmTMB::nbinom2,
#       priors = data.frame(
#         prior = prior,
#         class = class,
#         coef = coef) ))
# 
#   m1.reg_results <- seqwrap(
#     m1.reg,
#     return_models = FALSE,
#     verbose = FALSE,
#     #  subset = 1:50, # for prototyping
#     cores = parallel::detectCores())
# 
# 

#| label: load-sim-results
#| echo: false
#| message: false
#| warning: false





library(cowplot)
source("figures/figure-opts.R")


## Load data from simulations
source("R/simulation-functions.R")

# Extract simulations sorts out all simulations files and combines them.
# This is simulations with dispersion scenario 1.
sim_results <- extract_simulations(evaluations_path = "data_sim/evaluations",
                                   estimates_path = "data_sim/estimates",
                                   populationeffects_path = "data_sim/simdata/popeffect",
                                   disp_scenario = "s1")


sim_results2 <- extract_simulations(evaluations_path = "data_sim/evaluations2",
                                   estimates_path = "data_sim/estimates2",
                                   populationeffects_path = "data_sim/simdata2/popeffect",
                                   disp_scenario = "s2")



#| label: fig-genecounts-simulated
#| echo: false
#| message: false
#| warning: false
#| fig-cap: Number of simulated gene targets (true negative and true positives) across
#|   dispersion scenarios and sample sizes.
#| fig-height: 6



# Number of genes in each simulated data set
bind_rows(sim_results$populationeffects, 
          sim_results2$populationeffects) |> 
  # Filter to keep only 1 term
  filter(term == "conditionB") |> 
  mutate(effect = if_else(population_effect == 0, "true.negative", "true.positive")) |> 
  summarise(.by = c(dataset, size, disp_scenario, effect),
            n = n()) |>
  
  mutate(size = factor(size, levels = c("small", "medium", "large"), 
                       labels = c("Small (m = 54, n = 18)", 
                                  "Medium (m = 108, n = 36)", 
                                  "Large (m = 216, n = 72)")), 
         disp_scenario = factor(disp_scenario, 
                                levels = c("s1", "s2"), 
                                labels = c("High dispersion\nvariability", 
                                           "Low dispersion\nvariability")), 
         effect = factor(effect, levels = c("true.negative", "true.positive"), 
                         labels = c("True negative", "True positive"))) |> 
  
  ggplot(aes(size, n, fill = disp_scenario)) + 
  geom_point(position = position_jitter(width = 0.05), 
             size = 3, 
             shape = 21,
             alpha = 0.5) +
  labs(x = "Sample size (m samples, n participants)", 
       y = "Gene targets", 
       fill = "Dissersion scenario")  +
  theme_classic() + 
  
  
  
  scale_fill_manual(values = c(colors[c(1,5)])) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0)) + 
  facet_wrap(~ effect, scales = "free", ncol = 1)


#| echo: false
#| message: false
#| warning: false


# Number of genes per data set with true population effect
true_effects <- bind_rows(sim_results$populationeffects, 
          sim_results2$populationeffects) |> 
  mutate(effect = if_else(population_effect == 0, "true.negative", "true.positive")) |> 
  summarise(.by = c(term, dataset, size, disp_scenario, effect), 
            n = n()) 

# Savining for manuscript

# Check if data-out is a folder
if(!dir.exists("data-out")) dir.create("data-out")

saveRDS(true_effects, "data-out/true_effects.RDS")

# Combine all true/false effects
est_temp <- bind_rows(sim_results$estimates, 
          sim_results2$estimates) |> 
  # Retain only terms of interest.
  filter(term %in% c("conditionB", "timet3:conditionB")) |>
  # dataset was mis-named in the simulations.
  dplyr::rename(dataset = datasets) |>
  dplyr::select(- file)


est <- est_temp |> 
  inner_join(
    bind_rows(sim_results$populationeffects, 
              sim_results2$populationeffects) |>
               mutate(target = as.character(target)) |>
               dplyr::select(-file)) |>

  dplyr::select(target, term, estimate:p.value, population_effect, dataset, model, size, disp_scenario) |>

  mutate(.by = c(model, term, size, dataset, disp_scenario),
         fdr = p.adjust(p.value, method = "fdr")) |>
  mutate(true_effect = if_else(population_effect == 0, "neg", "pos"),
         identified_effect = if_else(fdr > 0.05, "neg", "pos"),
         true_positive = if_else(true_effect == "pos" &
                                   identified_effect == "pos", TRUE, FALSE),
         false_positive = if_else(true_effect == "neg" &
                                    identified_effect == "pos", TRUE, FALSE)) 



#| echo: false
#| message: false
#| warning: false
#| eval: true


convergence_stats <- bind_rows(sim_results$evaluations, 
          sim_results2$evaluations) |> 
  dplyr::select(-file) |> 
  mutate(conv = if_else(model == "m3" & isSingular == TRUE, FALSE, 
                        if_else(model != "m3" & convergence != 0, FALSE, TRUE)), 
         conv.nonstrict = if_else(model != "m3" & convergence != 0, FALSE, TRUE)) |> 
  dplyr::rename(dataset = datasets) |> 
   dplyr::select(target, size, model, dataset, disp_scenario, isSingular, conv, conv.nonstrict) |> 

  inner_join(
    bind_rows(sim_results$populationeffects, 
              sim_results2$populationeffects) |>
               mutate(target = as.character(target)) |>
                dplyr::select(-file) |> 
      mutate(effect = if_else(population_effect == 0, "n", "p")) |> 
       dplyr::select(target, term, dataset, size, disp_scenario, effect) |> 
      pivot_wider(names_from = term, values_from = effect) |> 
      mutate(effect = paste0(conditionB, ":", `timet3:conditionB`)) |> 
      mutate(.by = c(dataset, size, effect, disp_scenario), 
             ntotal = n()) |> 
       dplyr::select(-conditionB, - `timet3:conditionB`) 
  )  |> 
  
  summarise(.by = c(model, dataset, disp_scenario, effect, size), 
            conv = sum(conv),
            conv.nonstrict = sum(conv.nonstrict),
            ntotal = mean(ntotal)) 
# Save for manuscript  
saveRDS(convergence_stats, "data-out/covergence_stats.RDS")  


#| label: fig-convergence
#| echo: false
#| message: false
#| warning: false
#| eval: true
#| fig-cap: Number of models that provided estimates as a percentage of all potential
#|   targets. In (A), a strict criteria was used for convergence for the Gaussian model
#|   of transformed counts where models with singular fits where excluded as suggested
#|   by Vestal et al. (2022). In (B) all converged models are included.
#| fig-height: 8



p1 <- convergence_stats |> 

  mutate(success = 100 *  (conv / ntotal), 
         size = factor(size, levels = c("small", "medium", "large"), 
                       labels = c("Small\n(m = 54, n = 18)", 
                                  "Medium\n(m = 108, n = 36)", 
                                  "Large\n(m = 216, n = 72)")), 
         disp_scenario = factor(disp_scenario, 
                                levels = c("s1", "s2"), 
                                labels = c("High dispersion\nvariability", 
                                           "Low dispersion\nvariability")), 
         model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"), 
                        labels = c("Negative binomial", 
                                   "Regularized Negative binomial", 
                                   "Gaussian transformed counts", 
                                   "Poisson OLRE", 
                                   "Regularized Poisson OLRE"))) |> 
  ggplot(aes(model, success, color = model)) +
  geom_point(position = position_jitter(width = 0.1)) + 
  facet_grid(disp_scenario ~ size) +
  
  labs(y = "Converged/non-singular models (%)", 
       x = "Model") +
  
  theme_classic() + 
  theme(panel.background = element_rect(fill = "gray95"), 
        strip.background = element_blank(), 
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(hjust = 0), 
        legend.title = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = "bottom", 
        legend.box = "vertical") + 
  scale_color_manual(values = colors) + 
  guides(color=guide_legend(nrow=3, byrow=TRUE))
  
  
  
p2 <-  convergence_stats |> 

  mutate(success = 100 *  (conv.nonstrict / ntotal), 
         size = factor(size, levels = c("small", "medium", "large"), 
                       labels = c("Small\n(m = 54, n = 18)", 
                                  "Medium\n(m = 108, n = 36)", 
                                  "Large\n(m = 216, n = 72)")), 
         disp_scenario = factor(disp_scenario, 
                                levels = c("s1", "s2"), 
                                labels = c("High dispersion\nvariability", 
                                           "Low dispersion\nvariability")), 
         model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"), 
                        labels = c("Negative binomial", 
                                   "Regularized Negative binomial", 
                                   "Gaussian transformed counts", 
                                   "Poisson OLRE", 
                                   "Regularized Poisson OLRE"))) |> 
  ggplot(aes(model, success, color = model)) +
  geom_point(position = position_jitter(width = 0.1)) + 
  facet_grid(disp_scenario ~ size) +
  
  labs(y = "Converged models (%)", 
       x = "Model") +
  
  theme_classic() + 
  theme(panel.background = element_rect(fill = "gray95"), 
        strip.background = element_blank(), 
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(hjust = 0), 
        legend.title = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = "none")  + 
  scale_color_manual(values = colors) 
  
plot_grid(p1, 
          p2, 
          ncol = 1, 
          rel_heights = c(1, 0.8),
          labels = c("A", "B"))  


#| label: fig-pdist
#| echo: false
#| message: false
#| warning: false
#| fig-cap: Distribution of unadjusted p-values from all null-effects in models of simulated
#|   data. Non-uniform distributions indicate mis-specified models, i.e. models that
#|   do not capture variation in the data. Models are the negative binomial (NB), regularized
#|   negative binomial (NB-R), the Gaussian model of transformed counts (Gaussian), the
#|   Poisson observation-level random effect model (POLRE) and the regularized Poisson
#|   observation-level random effect model (POLRE-R).
#| fig-height: 8
#| fig-width: 6


# P-value distributions for null genes

est |>
  filter(population_effect == 0) |>
  
  mutate(model_coef = factor(paste0(model,"_", term), 
                             levels = c("m1_conditionB", 
                                        "m2_conditionB", 
                                        "m3_conditionB", 
                                        "m4_conditionB", 
                                        "m5_conditionB",
                                        "m1_timet3:conditionB", 
                                        "m2_timet3:conditionB",
                                        "m3_timet3:conditionB",
                                        "m4_timet3:conditionB",
                                        "m5_timet3:conditionB"
                                        ), 
                             labels = c("NB\nMain", 
                                        "NB-R\nMain",
                                        "Gaussian\nMain",
                                        "POLRE\nMain",
                                        "POLRE-R\nMain",
                                        "NB\nInteraction",
                                        "NB-R\nInteraction",
                                        "Gaussian\nInteraction",
                                        "POLRE\nInteraction",
                                        "POLRE-R\nInteraction"
                                        )), 
         size_disp = factor(paste0(size,"_", disp_scenario), 
                            levels = c("small_s2", "medium_s2", "large_s2", 
                                       "small_s1", "medium_s1", "large_s1"), 
                            labels = c("Small\nLow dispersion", 
                                       "Medium\nLow dispersion", 
                                       "Large\nLow dispersion", 
                                       "Small\nHigh dispersion", 
                                       "Medium\nHigh dispersion", 
                                       "Large\nHigh dispersion"
                                       ))) |> 
  
  ggplot(aes(p.value)) +
  geom_histogram(color="black", fill=colors[2], binwidth = 0.025,
                 boundary = 0, closed = "left", 
                 linewidth = 0.1) +
  scale_x_continuous(breaks = c(0.05, 0.5, 1)) +
  facet_grid(model_coef ~ size_disp) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 8), 
        axis.text = element_text(size = 8)) +
  labs(x = "Un-adjusted P-value", 
       y = "")




#| label: fig-cor-sim
#| echo: false
#| message: false
#| warning: false
#| fig-cap: Average correlations between simulated population effects and estimates as
#|   a function of sample size in two parameters (Main and interaction effects) and two
#|   dispersion scenarios.
#| fig-height: 5
#| fig-width: 6


# Model averages
mod_avg <- est |>

  filter(population_effect != 0) |>

  mutate(.by = c(model, size, term, dataset, disp_scenario),
         rank_true = rank(population_effect),
         rank_obs = rank(estimate)) |>

  summarise(.by = c(model, size, term, disp_scenario),
            cor = cor(rank_true, rank_obs)) |>

 mutate(size = factor(size, levels = c("small", "medium", "large"), 
                       labels = c("Small", "Medium", "Large")),
         model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                        labels = c("Negative binomial (non-informed)",
                                   "Regularized Negative binomial",
                                   "Gaussian transformed counts",
                                   "Poisson OLRE (non-informed)",
                                   "Regularized Poisson OLRE")), 
         disp_scenario = factor(disp_scenario, levels = c("s2", "s1"), 
                                labels = c("Low variability", 
                                           "High variability")), 
        term = factor(term, levels = c("conditionB", "timet3:conditionB"), 
                      labels = c("Main effect", 
                                 "Interaction effect"))) 



est |>

  filter(population_effect != 0) |>

  mutate(.by = c(model, size, term, dataset, disp_scenario),
         rank_true = rank(population_effect),
         rank_obs = rank(estimate)) |>

  summarise(.by = c(model, size, term, dataset, disp_scenario),
            cor = cor(rank_true, rank_obs)) |>

 mutate(size = factor(size, levels = c("small", "medium", "large"), 
                       labels = c("Small", "Medium", "Large")),
         model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                        labels = c("Negative binomial (non-informed)",
                                   "Regularized Negative binomial",
                                   "Gaussian transformed counts",
                                   "Poisson OLRE (non-informed)",
                                   "Regularized Poisson OLRE")), 
         disp_scenario = factor(disp_scenario, levels = c("s2", "s1"), 
                                labels = c("Low variability", 
                                           "High variability")), 
        term = factor(term, levels = c("conditionB", "timet3:conditionB"), 
                      labels = c("Main effect", 
                                 "Interaction effect")))  |>
  
  ggplot(aes(size, cor, group = paste(dataset, model),
             color = model)) +
  
    
  geom_line(alpha = 0.3) +
  
  geom_point(data = mod_avg, 
             color = "black",
             aes(group = NULL, shape = model, fill = model), 
             size = 3) +
  
  scale_fill_manual(values = c(colors[1], colors[2], colors[4], colors[1], colors[2])) +
  scale_color_manual(values = c(colors[1], colors[2], colors[4], colors[1], colors[2], 
                                colors[1], colors[2], colors[4], colors[1], colors[2])) +
  scale_shape_manual(values = c(22,22, 23, 25, 25)) +
  

  facet_grid(term ~ disp_scenario) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        legend.title = element_blank(), 
        axis.title.x = element_blank()) +
  labs(y = "Average correlation (r)")




#| echo: false
#| eval: false
# 
# est |>
# 
#   filter(population_effect != 0) |>
#   filter(dataset == 2, size == "medium") |>
# 
#   mutate(.by = c(model, term, dataset, disp_scenario),
#          rank_true = rank(population_effect),
#          rank_obs = rank(estimate)) |>
# 
#   ggplot(aes(rank_true, rank_obs)) +
#   geom_point(alpha = 0.3) +
#   facet_grid(paste(model) ~ paste(term, disp_scenario))
# 
# 
# 
# 

#| label: tbl-gsea
#| echo: false
#| message: false
#| warning: false
#| tbl-cap: Gene set enrichment analysis of fold-changes in all models. For each model
#|   the top (and bottom) ranked gene sets were selected based on -log~10~(FDR) &times;
#|   NES. The sign of the normalized effect size (NES) indicate enrichment at top (positive
#|   fold-change in T2D vs. control) and bottom (negative fold change) ranked genes.



library(gt)
library(stringr)


# To get the gsea results we need to run the figure 4 script
if(!file.exists("data-out/gsea_results.RDS")) source("figures/figure4.R")

gsea_results <- readRDS("data-out/gsea_results.RDS")


gsea_up <- gsea_results |> 
  
  slice_max(order_by = (-log10(p.adjust) * NES), n = 2, by = model) |> 
  dplyr::select(-core_enrichment) |> 
    mutate(model_descr = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                        labels = c("Negative binomial (non-informed)",
                                   "Regularized Negative binomial",
                                   "Gaussian transformed counts",
                                   "Poisson OLRE (non-informed)",
                                   "Regularized Poisson OLRE"))) 

gsea_down <- gsea_results |> 
  
  slice_max(order_by = -(-log10(p.adjust) * NES), n = 2, by = model) |> 
  dplyr::select(-core_enrichment) |> 
    mutate(model_descr = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                        labels = c("Negative binomial (non-informed)",
                                   "Regularized Negative binomial",
                                   "Gaussian transformed counts",
                                   "Poisson OLRE (non-informed)",
                                   "Regularized Poisson OLRE"))) 



bind_rows(gsea_up, gsea_down) |>  
          mutate(NES = round(NES, 2)) |>
  dplyr::select(model_descr,ID, Description, setSize, NES, p.adjust) |> 
  mutate(Description = str_to_title(Description)) |> 
  
 gt(groupname_col = "model_descr") |> 
  fmt_scientific(columns = p.adjust) |> 
  opt_table_font(size = 10) |> 
  cols_label("ID"=  "Gene ontology ID", 
             "setSize" = "Set size", 
             "p.adjust" = "FDR")








# Example of elapsed time using different modeling strategies. 
library(seqwrap)

set.seed(1) 
dat <- simcounts(n_genes = 100)

# Store data sets
d <- dat$data
md <- dat$metadata


mod1 <- seqwrap_compose(data = d, 
                metadata = md, 
                modelfun = stats::lm, 
                arguments = list(formula = y ~ x), 
                samplename = "sample")

mod2 <- seqwrap_compose(data = d, 
                metadata = md, 
                modelfun = lme4::lmer, 
                arguments = list(formula = y ~ x + (1|cluster)), 
                samplename = "sample")

mod3 <- seqwrap_compose(data = d, 
                metadata = md, 
                modelfun = glmmTMB::glmmTMB, 
                arguments = list(formula = y ~ x + (1|cluster)), 
                samplename = "sample")

results1 <- seqwrap(mod1, verbose = FALSE)
results2 <- seqwrap(mod2, verbose = FALSE)
results3 <- seqwrap(mod3, verbose = FALSE)

results1@elapsed_time
results2@elapsed_time
results3@elapsed_time


#| label: fig-timing
#| echo: false
#| message: false
#| warning: false
#| fig-cap: 'Comparing elapsed time between three different algorithms using the same
#|   basic model (negative binomial on simulate count data with 8 clusters and 2 observations
#|   each). '


dat <- simcounts(n_genes = 1000)

d <- dat$data
md <- dat$metadata

md$libsize <- log(colSums(d[,-1]))

available_cores <- parallel::detectCores()

size <- c(62, 125, 250, 500)
cores <- c(2, available_cores, available_cores, available_cores)

time_results <- data.frame(size = size, 
                      mod1 = rep(NA, length(size)), 
                      mod2 = rep(NA, length(size)), 
                      mod3 = rep(NA, length(size)))

for(i in seq_along(size)) {
  
  d_sub <- d[1:size[i], ]
  
  
  
  mod1 <- seqwrap_compose(data = d_sub, 
                metadata = md, 
                modelfun = lme4::glmer.nb, 
                arguments = list(formula = y ~ x + (1|cluster) + offset(libsize)), 
                samplename = "sample")
  
  mod2 <- seqwrap_compose(data = d_sub, 
                metadata = md, 
                modelfun = glmmTMB::glmmTMB, 
                arguments = list(formula = y ~ x + (1|cluster) + offset(libsize), 
                                 family = glmmTMB::nbinom2), 
                samplename = "sample")
  
  mod3 <- seqwrap_compose(data = d_sub, 
                metadata = md, 
                modelfun = glmmTMB::glmmTMB, 
                arguments = list(formula = y ~ x + (1|cluster) + offset(libsize), 
                                 family = glmmTMB::nbinom2, 
                                 control=glmmTMB::glmmTMBControl(optimizer=optim,
                                                        optArgs=list(method="BFGS"))), 
                samplename = "sample")
  
  
  
  
  
  
  results1 <- seqwrap(mod1, verbose = FALSE, cores = cores[i])
  results2 <- seqwrap(mod2, verbose = FALSE, cores = cores[i])
  results3 <- seqwrap(mod3, verbose = FALSE, cores = cores[i])
  
  
  time_results[i,2] <- results1@elapsed_time[3]
  time_results[i,3] <- results2@elapsed_time[3]
  time_results[i,4] <- results3@elapsed_time[3]

}



time_results |> 
  pivot_longer(cols = c(mod1, mod2, mod3)) |> 
  
  mutate(name = factor(name, levels = c("mod1", "mod2", "mod3"), 
                       labels = c("lme4::glmer.nb", 
                                  "glmmTMB (defualt optimizer)", 
                                  "glmmTMB optimizer: 'optim' - method 'BFGS'"))) |> 
  
  ggplot(aes(size, value, group = name, color = name)) + 
  geom_line() +
  geom_point() + 
  theme_classic() + 
  labs(color = "", 
       x = "Number of targets in data set", 
       y = "Computation time (sec)")
  
  




