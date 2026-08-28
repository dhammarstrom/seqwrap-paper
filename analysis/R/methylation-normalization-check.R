# This script asks two questions about the case study that the full model run
# cannot answer on its own.
#
# (1) Does the normalization change the conclusions? @seaborne2018 states that
#     the arrays were normalized with SWAN (subset-quantile within-array
#     normalization, Maksimovic et al. 2012), a within-array method that
#     equalises the type I and type II probe distributions inside each array
#     and does nothing across arrays. The pipeline here
#     (analysis/R/methylation-case-study-dataprep.R) uses functional
#     normalization, which uses the control probes to remove between-array
#     technical variation. The two differ in exactly the place where a
#     genome-wide shift between time points would live, so the direction
#     imbalance reported in the case study has to be checked against the
#     normalization it was computed on.
#
# (2) Do the positions the paper singles out behave the same way in our fits?
#     Figure 5 of @seaborne2018 names RPL35a, UBR5, SETD3, PLA2G16 and HEG1,
#     and the clustering in the same paper names a further eleven genes. No
#     probe identifiers are given in the text, so the genes are resolved to
#     probes through the array annotation, and every probe annotated to one of
#     them is reported rather than a selected one.
#
# The model fits are on a SUBSET: the mixed models cost hours over the whole
# array (see the timings in the case study), and a random subset of the size
# used here settles a question about the bulk behaviour of the arrays to well
# inside the precision anything is read at. The cheap tests -- the two-group
# comparison the published analysis used, and the genome-wide summaries of the
# beta values themselves -- are run over every position that survives
# filtering under BOTH normalizations, since they cost seconds.

## Packages and settings ###################################################
library(minfi)
library(dplyr)
library(seqwrap)

set.seed(2)   # preprocessSWAN() samples probes internally

cores <- parallel::detectCores()

out_dir <- here::here("analysis/data/derived_data/normalization_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

der_dir <- here::here("analysis/data/derived_data")

n_subset <- 10000   # random positions carried into the model fits

# The genes named in the paper. Figure 5 names the first five; the rest are
# the cluster A / cluster B members named in the results text. Matched against
# the semicolon-separated UCSC_RefGene_Name field of the array annotation,
# case-insensitively, since the paper writes RPL35a and the annotation RPL35A.
fig5_genes <- c("RPL35A", "UBR5", "SETD3", "PLA2G16", "HEG1")
cluster_genes <- c("C12ORF50", "BICC1", "ZFP2", "ODF2",
                   "AXIN1", "GRIK2", "CAMK4", "TRAF1", "NR2F6", "RSU1",
                   "STAG1")
paper_genes <- c(fig5_genes, cluster_genes)


## Inputs ##################################################################

rgset <- readRDS(file.path(der_dir, "seaborne-rgset.RDS"))
metadata <- readRDS(file.path(der_dir, "seaborne-metadata.RDS"))
gset_fun <- readRDS(file.path(der_dir, "seaborne-gset-normalized.RDS"))

stopifnot(setequal(colnames(gset_fun), metadata$geo_accession))
rgset <- rgset[, colnames(gset_fun)]


## SWAN, filtered the same way ##############################################
#
# The filtering is lifted from the dataprep script deliberately: comparing a
# SWAN set to the functional-normalization set is only a comparison of
# normalizations if everything downstream of it is identical.

swan_file <- file.path(out_dir, "gset-swan.RDS")

if (!file.exists(swan_file)) {
  det_p <- minfi::detectionP(rgset)

  gset_swan <- minfi::mapToGenome(minfi::preprocessSWAN(rgset))

  det_p <- det_p[match(minfi::featureNames(gset_swan), rownames(det_p)), ]
  gset_swan <- gset_swan[rowSums(det_p < 0.01) == ncol(gset_swan), ]
  gset_swan <- minfi::dropLociWithSnps(gset_swan, snps = c("SBE", "CpG"),
                                       maf = 0)
  gset_swan <- maxprobes::dropXreactiveLoci(gset_swan)
  anno_swan <- minfi::getAnnotation(gset_swan)
  gset_swan <- gset_swan[!anno_swan$chr %in% c("chrX", "chrY"), ]

  saveRDS(gset_swan, swan_file)
  rm(det_p)
  gc()
} else {
  gset_swan <- readRDS(swan_file)
}

rm(rgset)
gc()

message("positions after filtering: funnorm ", nrow(gset_fun),
        ", swan ", nrow(gset_swan))


## Beta values on both normalizations ######################################

beta_fun <- minfi::getBeta(gset_fun, offset = 100)
beta_swan <- minfi::getBeta(gset_swan, offset = 100)
beta_swan <- beta_swan[, colnames(beta_fun)]

common <- intersect(rownames(beta_fun), rownames(beta_swan))
beta_fun <- beta_fun[common, ]
beta_swan <- beta_swan[common, ]

message("positions common to both: ", length(common))

anno <- minfi::getAnnotation(gset_fun)
anno <- anno[common, ]

# Probes annotated to the genes the paper singles out.
gene_hit <- vapply(strsplit(toupper(anno$UCSC_RefGene_Name), ";"),
                   function(g) any(g %in% paper_genes), logical(1))
gene_probes <- rownames(anno)[gene_hit]

gene_map <- data.frame(
  id = rownames(anno)[gene_hit],
  gene = vapply(strsplit(toupper(anno$UCSC_RefGene_Name[gene_hit]), ";"),
                function(g) paste(sort(unique(g[g %in% paper_genes])),
                                  collapse = "/"),
                character(1)),
  group = NA_character_,
  chr = anno$chr[gene_hit],
  pos = anno$pos[gene_hit],
  island = anno$Relation_to_Island[gene_hit],
  region = vapply(strsplit(anno$UCSC_RefGene_Group[gene_hit], ";"),
                  function(g) paste(unique(g), collapse = "/"), character(1)),
  stringsAsFactors = FALSE
)
gene_map$group <- ifelse(
  vapply(strsplit(gene_map$gene, "/"), function(g) any(g %in% fig5_genes),
         logical(1)), "figure 5", "cluster")

message("probes annotated to the ", length(paper_genes), " named genes: ",
        nrow(gene_map))

saveRDS(gene_map, file.path(out_dir, "gene-map.RDS"))


## Genome-wide, both normalizations, no model ##############################
#
# Two summaries that cost nothing and do not depend on the subset:
#
# `shift`  the distribution of the per-position mean difference from baseline.
#          The case study reads a genome-wide direction imbalance off the
#          model estimates; this is the same quantity in beta units, before
#          any model, computed under both normalizations.
#
# `twogrp` the published test itself -- a two-group comparison of baseline
#          against the later time point with participants treated as
#          independent -- under both normalizations, so that the count of
#          positions it returns can be attributed to the test or to the
#          normalization.

lev <- c(time3_loading = "3_loading", time4_unloading = "4_unloading",
         time5_reloading = "5_reloading")
tm <- metadata$time[match(colnames(beta_fun), metadata$geo_accession)]
i_base <- which(tm == "1_baseline")

two_group <- function(mat, i1, i2) {
  n1 <- length(i1); n2 <- length(i2)
  m1 <- rowMeans(mat[, i1, drop = FALSE])
  m2 <- rowMeans(mat[, i2, drop = FALSE])
  ss <- rowSums((mat[, i1, drop = FALSE] - m1)^2) +
    rowSums((mat[, i2, drop = FALSE] - m2)^2)
  df2 <- n1 + n2 - 2
  tstat <- (m1 - m2) / sqrt((ss / df2) * (1 / n1 + 1 / n2))
  data.frame(id = rownames(mat), diff = m2 - m1,  # condition minus baseline
             p = 2 * pt(-abs(tstat), df2))
}

norm_sets <- list(funnorm = beta_fun, swan = beta_swan)

twogrp <- bind_rows(lapply(names(norm_sets), function(nm)
  bind_rows(lapply(names(lev), function(t) {
    r <- two_group(norm_sets[[nm]], i_base, which(tm == lev[[t]]))
    r$term <- t
    r$norm <- nm
    r
  }))))

saveRDS(twogrp, file.path(out_dir, "twogroup.RDS"))

shift <- twogrp |>
  summarise(.by = c(norm, term),
            n = n(),
            median_diff = median(diff),
            mean_diff = mean(diff),
            prop_hypo = mean(diff < 0),
            p05 = sum(p < 0.05),
            prop_hypo_p05 = mean(diff[p < 0.05] < 0),
            fdr05 = sum(p.adjust(p, "BH") < 0.05))

saveRDS(shift, file.path(out_dir, "shift.RDS"))
print(as.data.frame(shift))


## The model subset ########################################################

subset_ids <- union(sample(common, min(n_subset, length(common))),
                    gene_probes)
message("positions carried into the model fits: ", length(subset_ids))

mk_df <- function(mat, ids) {
  m <- log2(mat[ids, ] / (1 - mat[ids, ]))   # M-values
  data.frame(id = rownames(m), m, check.names = FALSE)
}

fit_dfs <- list(funnorm = mk_df(beta_fun, subset_ids),
                swan = mk_df(beta_swan, subset_ids))

rm(beta_fun, beta_swan, gset_swan, gset_fun)
gc()


## Fits ####################################################################
#
# The m_reml arm of the full run, unchanged: gaussian likelihood on M-values,
# REML, read with both Wald reference distributions. It is the arm the
# permutation study found to hold its nominal size, so it is the one a
# conclusion should be read on. eval_fun and summary_fun are copied from
# analysis/R/methylation-case-full.R rather than re-derived.

eval_fun <- function(m) {
  vc <- glmmTMB::VarCorr(m)$cond$participant
  data.frame(convergence = m$fit$convergence,
             singular = !m$sdr$pdHess,
             re_sd = sqrt(as.numeric(vc[1, 1])),
             n = nrow(m$frame),
             reml = isTRUE(m$modelInfo$REML))
}

summary_fun_wald <- function(m) {
  cf <- coef(summary(m))$cond
  out <- data.frame(term = rownames(cf), estimate = cf[, "Estimate"],
                    std.error = cf[, "Std. Error"], p_wald = cf[, "Pr(>|z|)"],
                    row.names = NULL, stringsAsFactors = FALSE)
  sat <- tryCatch(coef(summary(m, ddf = "satterthwaite"))$cond,
                  error = function(e) NULL)
  i <- if (is.null(sat)) NA_integer_ else match(out$term, rownames(sat))
  out$ddf <- if (is.null(sat)) NA_real_ else sat[i, "ddf"]
  out$p_satt <- if (is.null(sat)) NA_real_ else sat[i, "Pr(>|t|)"]
  out
}

for (nm in names(fit_dfs)) {
  f <- file.path(out_dir, sprintf("fit-%s.RDS", nm))
  if (file.exists(f)) {
    message("skipping ", nm, " -- ", basename(f), " exists")
    next
  }

  message("fitting ", nm, " on ", cores, " cores")

  container <- seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    data = fit_dfs[[nm]],
    metadata = metadata,
    samplename = "geo_accession",
    eval_fun = eval_fun,
    summary_fun = summary_fun_wald,
    arguments = list(formula = y ~ time + (1 | participant), REML = TRUE))

  results <- seqwrap(container, return_models = FALSE, cores = cores,
                     verbose = FALSE)
  summarised <- seqwrap_summarise(results, verbose = FALSE)

  saveRDS(list(norm = nm,
               summaries = summarised$summaries,
               evaluations = summarised$evaluations,
               errors = summarised$errors,
               elapsed_time = results@elapsed_time,
               k = results@k, n = results@n),
          f)

  message(nm, " done in ", round(results@elapsed_time[["elapsed"]] / 60, 1),
          " minutes")

  rm(container, results, summarised)
  gc()
}


## Collected model results #################################################

fits <- bind_rows(lapply(names(fit_dfs), function(nm) {
  x <- readRDS(file.path(out_dir, sprintf("fit-%s.RDS", nm)))
  x$summaries |>
    inner_join(select(x$evaluations, target, convergence, singular, re_sd),
               by = "target") |>
    mutate(norm = nm)
}))

saveRDS(fits, file.path(out_dir, "fits.RDS"))

message("model fits: ", nrow(fits), " rows over ",
        n_distinct(fits$target), " positions")
