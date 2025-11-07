library(glmmSeq)
library(dplyr)
library(tidyverse)
library(edgeR)
library(seqwraphelper)



dat <- geo_data()

#
data_df <- dat$countdata %>%
  # filter rows with a zero variance
  filter(apply(dplyr::select(., -1), 1, var) != 0) %>%
  distinct(geneid, .keep_all = T)

#saveRDS(gene_data, "data/alcohol_study_gene_data.RDS")

metadata <- dat$metadata %>%
  # subset to include only those whose sequence data are available
  filter(seq_sample_id %in% colnames( data_df[, -1])) %>%
  # drop the tissue colume
  dplyr::select(-tissue) %>%
  mutate(sex = factor(sex, levels =c("Male", "Female")),
         dose = factor(dose, levels = c("Placebo", "MediumDose", "HighDose")),
         period = factor(period, levels = c("1", "2","3")),
         time = factor(time, levels = c("0", "3")),
         study_site = factor(study_site, levels = c("1", "2")))


#saveRDS(metadata, "data/alcohol_study_metadata.RDS")

# Check if everything matches except the gene_id
all(colnames( data_df)[-1] == metadata$seq_sample_id)

 # we would be using the steps shown on glmmseq's vignette https://myles-lewis.github.io/glmmSeq/articles/glmmSeq.html

 # estimate dispersion using edgeR

 disp <- setNames(edgeR::estimateDisp(data_df[,-1])$tagwise.dispersion, rownames(data_df))

 # Monitor run time
# start_time <- Sys.time()


 glmmseq_model <- glmmSeq( ~ dose*period +time + sex + study_site + (1 | participant),
                          countdata = data_df[,-1],
                          metadata = metadata,
                          dispersion = disp,
                          cores = 4)


 # end_time <- Sys.time()
 # elapsed <- end_time - start_time
 # print(elapsed)



 formula <- age ~ dose + gender
 glm(formula, data = metadata)

