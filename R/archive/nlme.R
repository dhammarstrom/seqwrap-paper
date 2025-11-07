# This script uses SeqWrap to build a model using nlme
library(seqwrap)
library(dplyr)
library(tidyverse)
library(nlme)
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
length(unique(metadata$participant))
# Check if everything matches except the gene_id
all(colnames(gene_data)[-1] == metadata$seq_sample_id)



seqwrap:::

# Check if everything matches except the transcript_id
match(colnames(gene_data), metadata$seq_sample_id)

# test and see how the model runs in lme
dx <- gene_data %>%
  slice_head(n= 1) %>%
  pivot_longer(cols = -geneid, names_to = "seq_sample_id", values_to = "count") %>%
  inner_join(metadata, by = "seq_sample_id")


fit_nlm <- lme(fixed = count ~ dose + period + time + sex + study_site  ,
               random = ~ 1 | participant,
               data = dx)

x <- summary(fit_nlm)

x$corFixed # shows the correlation between the fixed effects
x$tTable # fixed effects table
x$modelStruct$reStruct # rnadom effects
x$fitted #  fitted values represent the model’s predicted values for the response variable,
#based on the fixed effects and the estimated random effects

sum_fun_lme <- function(x){

  data.frame(coef(summary(x)))

}

seqwrap:::simcounts() #simulate counts for a specific type
# suggest workflow to save the models at first to check if it works

# additional vars

nlme_container <- seqwrap_compose(data = gene_data,
                             metadata = metadata,
                             samplename = "seq_sample_id",
                             modelfun = nlme::lme ,
                             arguments = list(fixed = y ~ dose + period + time + sex + study_site  ,
                                              random = ~1|participant, na.action = na.omit),
                            # summary_fun = sum_fun_lme,
                             #additional_vars = metadata$age,
                           # eval_fun = evalfun
                            )


nlme_model <- seqwrap(nlme_container,

                 return_models = FALSE,
                 subset = 2:15,
                 cores = 1)





model_summary <- seqwrap_summarise(nlme_model)

model_summary$summaries
