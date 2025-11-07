library(lmerseq)
library(dplyr)
library(tidyverse)
library(edgeR)
library(seqwraphelper)

devtools::install_github("stop-pre16/lmerSeq", build_vignettes = TRUE)
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
