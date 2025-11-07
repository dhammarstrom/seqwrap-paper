library(seqwrap)
library(dplyr)
library(tidyverse)
library(lme4)
library(seqwraphelper)




# Load the gene data
gene_data <- readRDS("data/alcohol_study_gene_data.RDS")
# Load metadata
metadata <- readRDS("data/alcohol_study_metadata.RDS")



# Check if everything matches except the transcript_id
match(colnames(gene_data), metadata$seq_sample_id)

# test and see how the model runs in lme
dx <- gene_data %>%
  slice_head(n= 2) %>%
  pivot_longer(cols = -geneid, names_to = "seq_sample_id", values_to = "count") %>%
  inner_join(metadata, by = "seq_sample_id")


fit_lm <- glmer(count ~ dose * period * time + sex + study_site + (1 | participant),
               family = poisson, data = dx)

x <- summary(fit_lm)


# for coninous data, replace glmer with lmer
citation("lmerTest")

lmer_container <- seqwrap_compose(
  data = gene_data,
  metadata = metadata,
  samplename = "seq_sample_id",
  modelfun = lme4::glmer,
  arguments = list(
    formula = y ~ dose*time + sex +  period + (1+ time || participant) + (1|study_site),
  family = poisson),
#  summary_fun = sum_fun_lmer1,
 # eval_fun = evalfun
)

lmer_model <- seqwrap(
  lmer_container,
#  summary_fun = sum_fun_lmer1,
 # eval_fun = evalfun,
  return_models = F,
  subset = 1:15,
  cores = 2
)
x <- seqwrap_summarise(lmer_model)
x$summaries
saveRDS(lmer_model, "data/seqwrap_lmer_model.RDS")
