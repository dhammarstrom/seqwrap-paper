library(dplyr)
library(tidyverse)
library(stringr)

#read tsv file from the alcohol study https://pmc.ncbi.nlm.nih.gov/articles/PMC11087488/
# data downlaoded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE232408


# This contains the raw counts and
df <- readr::read_tsv("data/GSE232408_raw_counts_GRCh38.p13_NCBI.tsv")

# Load the the metadata file
lines <- readLines("data/GSE232408_series_matrix.txt")

# Filter lines starting with '!'
meta_lines <- lines[grepl("^!", lines)]

#Split each line into a list of values
meta_list <- lapply(meta_lines, function(line) {
  strsplit(line, "\t")[[1]]
})

# Extract keys and values
keys <- sapply(meta_list, function(x) x[1])
values <- lapply(meta_list, function(x) x[-1])

# Find the maximum length and pad shorter vectors
max_len <- max(sapply(values, length))
values_padded <- lapply(values, function(x) {
  length(x) <- max_len
  return(x)
})

# Convert to a data frame
metadata_df <- as.data.frame(do.call(cbind, values_padded), stringsAsFactors = FALSE)
colnames(metadata_df) <- gsub("^!", "", keys)
# there are duplicated column names. Make each unique
colnames(metadata_df) <- make.names(colnames(metadata_df), unique = T)


metadata_df <- metadata_df %>%
  # select the columns of interest based on index position
  dplyr::select(35,36, 46, 49, 50 ) %>%
  rename(participant_id = "Sample_title",
         sample_id = "Sample_geo_accession",
         dose = "Sample_characteristics_ch1.2",
         gender = "Sample_characteristics_ch1.5",
         age = "Sample_characteristics_ch1.6") %>%
  # extract the age
  mutate(age = as.numeric(str_extract(age, "\\d+\\.?\\d*"))) %>%
  # extract the experimental sequence from the participant_id
  separate(participant_id, into = c("participant_id", "experimental_sequence"), sep = "_") %>%
  separate(dose, into = c("x", "dose"), sep = ": ") %>%
  separate(gender, into = c("y", "gender"), sep = ": ")%>%
  select(-c(x,y))


# remove the quotation marks in the variable names

metadata_df[] <- lapply(metadata_df, function(x) {
  if (is.character(x)) gsub('"', '', x) else x
})

metadata_df <- metadata_df %>%
  mutate(gender = factor(gender, levels = c("Male", "Female")),
         dose = factor(dose, levels = c("Placebo", "MediumDose", "HighDose")),
         sequence = case_when(str_detect(experimental_sequence, "D0") ~ "D0",
                              str_detect(experimental_sequence, "D3") ~ "D3"),
         sequence = factor(sequence, levels = c("D0", "D3")))

length(unique(metadata_df$participant_id))

length(unique(metadata_df$experimental_sequence))

length(unique(metadata_df$sample_id))

length(unique(metadata_df$dose))


length(unique(metadata_df$age))

# There are only 53 samples in the gene expression data, and 70 in the metadata
# filter the metadata to include only samples in the gene expression

metadata <- metadata_df %>%
  filter(sample_id %in% colnames(df[, -1]))

# check if the order of sample names match in both dataframes

match(colnames(df), metadata$sample_id)

# match the gene_expression data's sample id to the metadata
# df <- df[, c("GeneID", metadata$sample_id)]

# filter rows with all zeros in the gene expression dataframe
nonZero <- df %>%
  filter(rowSums(df[, -1]) != 0)

# normaluise the gene expression data using TMM-normalised CPM values

# create dge list and calulate the normalisation factor
dge <- DGEList(nonZero[,-1])
dge <- calcNormFactors(dge, method = "TMM")

# get cpm values
data_df <- as.data.frame(cpm(dge, normalized.lib.sizes = T))
data_df$GeneID <- nonZero$GeneID
# make the GeneId the first column
data_df <- data_df %>%
  select(GeneID, everything()) %>%
  # round to 2 decimal places
  mutate(across(where(is.numeric), ~round(.x, 2)))





# save metadata and gene expression file

saveRDS(metadata, "data/alcohol_study_metadata.RDS")

saveRDS(data_df, "data/alcohol_study_gene_data.RDS")

