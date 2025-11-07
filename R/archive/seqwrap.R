library(seqwrap)
library(dplyr)
library(tidyverse)
library(stringr)
library(glmmSeq)
library(lme4)
library(nlme)
library(edgeR)
library(brms)


sumfun <- function(x) {
  # A function that returns a data frame with the coefficients
  # from the fitted model
  cond_effects <- data.frame(
    cbind(
      data.frame(
        coef = rownames(coef(summary(x))$cond)
      )
    ),
    coef(summary(x))$cond,
    row.names = NULL
  )

  return(cond_effects)
}

evalfun <- function(x) {
  # A function that returns a data frame with the p-values from
  # the residual diagnostics
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  results <- data.frame(
    pval.unif = unif$p.value,
  )

  return(results)
}
evalfun2 <- function(x) {
  # A function that returns a data frame with the p-values from
  # the residual diagnostics
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  results <- data.frame(
    pval.unif = unif$p.value
  )

}


# load the downloaded normalised counts
gene_data <- readr::read_tsv("data/GSE232408_norm_counts_TPM_GRCh38.p13_NCBI.tsv") %>%
  mutate(across(-1, ~ round(.x, 2))) %>%
  filter(apply(dplyr::select(., -1), 1, var) != 0)


metadata <- readRDS("data/alcohol_study_metadata.RDS")
# gene annoation downloaded  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE232408
 # same annotation used by study

 annotation <- readr::read_tsv("data/Human.GRCh38.p13.annot.tsv") %>%
   # select the gene_id, genetype and biotype
   dplyr::select(GeneID, Symbol, GeneType)



 gene_data <- gene_data %>%
   inner_join(annotation, by = "GeneID") %>%
   dplyr::select(-GeneID) %>%
   dplyr::select(Symbol, GeneType, everything())%>%
   dplyr::select(-GeneType)

 saveRDS(gene_data, "data/alcohol_gene_data_TPM.RDS")

 unique(gene_data$GeneType)

# subset the protein coding genes

 prot_coding <- gene_data %>%
   filter(GeneType == "protein-coding") %>%
   dplyr::select(-GeneType)
# select the non protein coding genes
 non_prot_coding <- gene_data %>%
   filter(GeneType == "ncRNA")%>%
   dplyr::select(-GeneType)


 saveRDS(prot_coding, "data/protein_coding_genes.RDS")
saveRDS(non_prot_coding, "data/ncRNAs.RDS")


# Transpose the protein-coding data

df_prot <- as.data.frame(t(prot_coding[, -1]))
# get the column names
colnames(df_prot) <- prot_coding[[1]]
# get the participant ids
df_prot$participant_id <- rownames(df_prot)




# Repeat for ncRNA

df_ncRNA <- as.data.frame(t(non_prot_coding[, -1]))

colnames(df_ncRNA) <- non_prot_coding[[1]]
# get the participant ids
df_ncRNA$participant_id <- rownames(df_ncRNA)












# To use multivariate mixed modelling, the two dependent variables should be coverted to lists of the same length

# To simulate that, we use the raw counts versus the TPM gene counts

TPM_counts <-  readr::read_tsv("data/GSE232408_norm_counts_TPM_GRCh38.p13_NCBI.tsv") %>%
  mutate(across(-1, ~ round(.x, 2))) %>%
  # filter(apply(dplyr::select(., -1), 1, var) != 0)
  filter(rowSums(df[, -1]) != 0) %>%
  # subset to the first 10 rows
  slice_head(n= 15)

raw_counts <- readRDS("data/alcohol_study_gene_data.RDS")%>%
  filter(apply(dplyr::select(., -1), 1, var) != 0) %>%
  # subset to the first 10 rows
  slice_head(n= 15)

# convert the two dataframes into a list
df_data <- list(TPM = TPM_counts, raw = raw_counts)
#names(df_data) <- c("protein_coding_genes", "ncRNA")


dx <- raw_counts %>%
  slice_head(n= 2) %>%
  pivot_longer(cols = -GeneID, names_to = "sample_id", values_to = "count") %>%
  inner_join(metadata, by = "sample_id")



fit_brms <- brm(count ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence),
                data = dx)


summary(fit_brms)


formula1 <- bf( age ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence))

formula2 <- bf (count ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence))

multi_fit_brms <- brm(formula1 + formula2 + set_rescor(T),
                      data = dx, family = gaussian())
fixef(multi_fit_brms)
multi_fit_brms$

x <- summary(multi_fit_brms)

x$fixed
x$formula
x$spec_pars
x$iter
x$algorithm
head(predict(multi_fit_brms))


plot(conditional_effects(multi_fit_brms))

pp_check(multi_fit_brms)



fit_nlm <- lme(count ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence),
               data = dx)


fit_lm <- lmer(count ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence),
             data = dx)

y <- summary(fit_lm)
y$coefficients



fit_glmm <- glmmTMB::glmmTMB(count ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence),
                 data = dx)

q <- summary(fit_glmm)

q$coefficients$cond
container <- seqwrap_compose(data = raw_counts,
                             metadata = metadata,
                             samplename = "sample_id",
                             modelfun = brms::brm,
                             arguments = list(formula = y ~ dose*sequence + gender + (1|participant_id) + (1|experimental_sequence),
                                              family = gaussian()),
                             summary_fun = sum_fun_brms,
                            # eval_fun = evalfun
                            )


model <- seqwrap(container,
                 summary_fun = sum_fun_brms,
                # eval_fun = evalfun2,
                 return_models = T,
                 subset = 1:150,
                 cores = 2)


model_summary <- seqwrap_summarise(model)

model
head(model_summary)
str(model)
model

view(seqwrap_compose)
view(seqwrap)
modelfun
attributes(container)
container@modelfun
