library(seqwrap)
library(tidyverse)
library(dplyr)


# This would generate a very basic control-study design

sim_study <- seqwrap:::simcounts(n_genes = 7000, n_samples = 100, clusters = 10)
sim_counts <- sim_study$data


sim_metadata <- sim_study$metadata %>%

  mutate( x = factor(x, levels = c("0", "1")))


# Check if everything matches except the transcript_id
match(colnames(sim_counts), sim_metadata$sample)

# create a container for an nlmelme model


sim_nlme_container <- seqwrap_compose(data = sim_counts,
                                  metadata = sim_metadata,
                                  samplename = "sample",
                                  modelfun = nlme::lme ,
                                  arguments = list(fixed = y ~ x ,
                                                   random = ~1|cluster,
                                                   na.action = na.omit))

# build the model
sim_nlme_model <- seqwrap(sim_nlme_container,

                      return_models = FALSE,
                      #subset = 2:15,
                      cores = 5)


model_summary <- seqwrap_summarise(sim_nlme_model)

nlme_sum <- model_summary$summaries

saveRDS(nlme_sum, "data/nlme_simulated_model_summary.RDS")





# simulation using  lme4


sim_lmer_container <- seqwrap_compose(data = sim_counts,
  metadata = sim_metadata,
  samplename = "sample",
  modelfun = lme4::glmer,
  arguments = list(
    formula = y ~ x + ( 1 | cluster),
    family = poisson),
  additional_vars = "cluster"
)



sim_lmer_model <- seqwrap(
  sim_lmer_container,
  return_models = F,
  #subset = 1:15,
  cores = 5)

lmer_model_summary <- seqwrap_summarise(sim_lmer_model)

lmer_sum <- lmer_model_summary$summaries

saveRDS(lmer_sum, "data/lmer_simulated_model_summary.RDS")



# simulation using glmmTMB

sim_glmm_container <- seqwrap_compose(data = sim_counts,
                                      metadata = sim_metadata,
                                      samplename = "sample",
                                      modelfun = glmmTMB::glmmTMB,
                                      arguments = list(
                                        formula = y ~ x + ( 1 | cluster),
                                        family = glmmTMB::nbinom1)
)

sim_glmm_model <- seqwrap(
  sim_glmm_container,
  return_models = F,
  #subset = 1:15,
  cores = 5)

glmm_model_summary <- seqwrap_summarise(sim_glmm_model)

glmm_sum <- glmm_model_summary$summaries

saveRDS(glmm_sum, "data/lmer_simulated_model_summary.RDS")
