# This script uses SeqWrap to build a model using glmmTMB

# Load needed libraries
library(seqwrap)
library(dplyr)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(seqwraphelper)
# This script uses lme in SeqWrap to analyse the alcohol -study data


# download test data using seqwraphelper
# this downloads a list containing the count data and the metadata
dat <- geo_data()

#
gene_data <- dat$countdata %>%
  # filter rows with a zero variance
  filter(apply(dplyr::select(., -1), 1, var) != 0) %>%
  distinct(geneid, .keep_all = T)

#saveRDS(gene_data, "data/alcohol_study_gene_data.RDS")

metadata <- dat$metadata %>%
  # subset to include only those whose sequence data are available
  filter(seq_sample_id %in% colnames(gene_data[, -1])) %>%
  # drop the tissue colume
  dplyr::select(-tissue) %>%
  mutate(sex = factor(sex, levels =c("Male", "Female")),
         dose = factor(dose, levels = c("Placebo", "MediumDose", "HighDose")),
         period = factor(period, levels = c("1", "2","3")),
         time = factor(time, levels = c("0", "3")),
         study_site = factor(study_site, levels = c("1", "2")))


#saveRDS(metadata, "data/alcohol_study_metadata.RDS")

# Check if everything matches except the gene_id
all(colnames(gene_data)[-1] == metadata$seq_sample_id)

# maeginal effects using mid to high dose


# initialise the argument
# negative binomial will be used
dose * time + period + study + (1| participant) + (1|site)

args <- list(formula = y ~ dose*period*time + sex + study_site + (1 | participant),
             family = glmmTMB::nbinom1)


glmm_container <- seqwrap_compose(data = gene_data,
                             metadata = metadata,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB ,
                             arguments = args,
                             summary_fun = sumfun,
                            eval_fun = evalfun)

glmm_model <- seqwrap(glmm_container,
                 summary_fun = sumfun,
                 eval_fun = evalfun,
                 return_models = F,
              #   subset = 1:1900,
                 cores = 2)




saveRDS(glmm_model, "data/seqwrap_glmm_model.RDS")
model_summary <- seqwrap_summarise(glmm_model)

model
head(model_summary)
str(glmm_model)
model
