
library(seqwrap)
library(dplyr)
library(tidyverse)
library(brms)
library(seqwraphelper)

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





# Compose the container with brms
# Define the wrapper call for brms
brms_container <- seqwrap_compose(
  data = gene_data,
  metadata = metadata,
  samplename = "seq_sample_id",
  modelfun = brms::brm,
  arguments = list(
    formula = y ~ dose * period * time + sex + study_site + (1 | participant),
    data = NULL,  # will be injected per gene
    family = brms::negbinomial(),
    chains = 4,
    cores = 4,
    iter = 2000,
    save_pars = brms::save_pars(all = TRUE),
    silent = TRUE
  )
  #,
 # summary_fun = sum_fun_brms,
#  eval_fun = evalfun_brms # assuming you have this defined
)


# Run the models
brms_model <- seqwrap(
  brms_container,
 # summary_fun = sum_fun_brms,
#  eval_fun = evalfun_brms,
  return_models = FALSE,
  subset = 1:15,
  cores = 2)


x <- seqwrap_summarise(brms_model)

x
