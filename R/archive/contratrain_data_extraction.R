

library(dplyr)
library(trainomeMetaData)
library(tidyverse)
library(biomaRt)
library(edgeR)
library(scales)



# Download metadata from  the TrainomeMetadata package
data("ct_participants")
data("ct_samples")  
data("ct_strength")
data("ct_thickness")


# rename the values in the time column to match those in the other dataframes
ct_strength["time"][ct_strength["time"] == "pre"] <- "t1"
ct_strength["time"][ct_strength["time"] == "mid"] <- "t3"
ct_strength["time"][ct_strength["time"] == "post"] <- "t4"

ct_metadata <- ct_samples %>% inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
  dplyr::select(study, participant, sex, condition,leg, time, seq_sample_id) %>%
  inner_join(ct_thickness, by = c("study", "participant", "leg", "condition", "time") ) %>%
  # drop the rows containing missing data
  drop_na()%>%
  # drop columns not needed for lncRNA analyses
  dplyr::select(-c("study", "leg", "thickness")) %>%
  # inner_join(ct_strength, by = c("study", "participant", "group", "leg", "condition", "time"))
  mutate(time = factor(time, levels = c("t1", "t3",   "t4")), 
         condition = factor(condition, levels = c("set0", "set3", "set6")), 
         sex = factor(sex, levels = c("male", "female")),
         training_status = if_else(condition == "set0", "untrained", "trained"),
         training_status = factor(training_status,
                                  levels = c("untrained", "trained"))
         ) 



# Visualize metadata


ct_metadata %>%
  ggplot(aes(condition, color = condition, fill = condition))+
  geom_bar( ) +
  facet_grid(time ~ group) +
  ggtitle("Distribution of RNA sequence data in the Contratrain study")+
  theme(axis.text = element_text(size = 10), text = element_text(size = 10),
        plot.title = element_text(hjust = 0.1))+
  ylab("Number of datapoints by time")+
  scale_y_continuous(labels = label_number(accuracy = 1))

#save figure in the plot folder
ggsave("plots/Data_distribution.png", width = 7, height = 7)





# download the contratrain gene counts from the trainomeMetadata
download_ome(download = "ct_gene_rsem")

# This would be used for model building
# It contains the expected counts of the genes
ct_genes <- read_csv("ome-data/ct_gene_rsem.csv") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  # drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))




# Load the RSEM FPKM counts fro the contratrain counts

download_ome(download = "ct_fpkm_rsem")

# This will be used during gene coexpression analyses
fpkm <- read_csv("ome-data/ct_fpkm_rsem.csv")%>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))


saveRDS(fpkm, file = "data/Ct_genes_fpkm.RDS")








# # Keep nonzero rows
nonzero <- ct_genes %>%
  dplyr::filter(rowSums(ct_genes[,-1]) != 0)
 
 
# # Keep filtered genes (based on group)
filtered <- nonzero %>%
  dplyr::filter(filterByExpr(nonzero[,-1], group = ct_metadata$time)) 


#saveRDS(filtered, file = "./data/filtered_all_gene_counts.RDS")




# Create dge lists and calculate norm factors
dge_ct   <- DGEList(filtered[,-1])


dge_ct  <- calcNormFactors(dge_ct, method = "TMM"  )


# Get the normalised counts

normalised_counts <- as.data.frame(cpm(dge_ct, normalized.lib.sizes = T))

normalised_counts$gene_name <- filtered$gene_name

# Make the gene_name the first column

normalised_counts <- normalised_counts %>%
  dplyr::select(gene_name, everything())

 df %>%
  filter(target %in% target_list)
print(subset_df)
#select the genes of interest from the dataframe
genes_of_int <- c("ENSG00000279762", "ENSG00000286733", "ENSG00000276564",
                  "PRKACB-DT", "MEF2C-AS1", "ATXN1-AS1",
                  "LYPLAL1-DT", "PKN2-AS1")

lncs_of_int <-  normalised_counts  %>%
  filter(gene_name %in% genes_of_int)%>%
  mutate_if(is.numeric, round, 2)

lncs <- lncs_of_int %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "count",
               cols = -(gene_name)) %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  group_by(gene_name, time) %>%
  summarise(average_count = mean(count))%>%
  mutate_if(is.numeric, round, 2) %>%
  pivot_wider(names_from = "time", values_from = "average_count")

writexl::write_xlsx(lncs, "data/DE_lncs.xlsx")





## Add effective library size to 
ct_metadata <- dge_ct$samples %>%
  rownames_to_column(var = "seq_sample_id") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  
  mutate(efflibsize = (lib.size * norm.factors) / median(lib.size * norm.factors)) %>%
  # Remove the group.x column as it was used for filtering
  dplyr::select(-("group.x")) %>%
  # rename the group.y column as it it the original grouping of the data
  # rename the library size and the normalization factors to reflect that all biotypes are accounted for
  rename("group" = "group.y") %>%
  
  print()

# Save metadata for further use

saveRDS(ct_metadata, file = "data/contratrain_metadata.RDS")



# Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

# extract GENE biotypes and  names
annotation<- getBM(attributes = c("gene_biotype", "external_gene_name"), mart = ensembl )
annotation <- inner_join(filtered, annotation, by= c("gene_name" = "external_gene_name"))%>%
  # Round the counts to zero decimal place
  mutate_at(2:137, ~ as.integer(., 0))


# check the unique biotypes
unique(annotation$gene_biotype)


# Extract only the lncRNAs and save in the data file


lncRNAS <- annotation %>%
  dplyr::filter(gene_biotype == "lncRNA") %>%
  dplyr::select(-gene_biotype)


# save the lncRNA count file
saveRDS(lncRNAS, file = "data/lncRNA_genes.RDS")



# Extract only the protein-coding genes
# Needed for correlation analyses

mRNAs <- annotation %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(-gene_biotype)


# Save the file
saveRDS(mRNAs, file = "data/protein_coding_genes.RDS")


# Using the already done annotation, select the fpkm counts of the mRNA genes 
prot <- fpkm[fpkm$gene_name %in% mRNAs$gene_name, ]

# Check for differences
setdiff(prot$gene_name, mRNAs$gene_name)

saveRDS(prot, "data/protein_coding_genes_FPKM.RDS")



