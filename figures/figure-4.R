# Gene set enrichment analysis ################################################
#
# 01. Load packages and data
# 02. Model performance (n estimates, AIC comp)
# 02. Filter genes with estimates per method (universe for enrichemnt analysis)
# and calculate statistics (FDR)
# 03. Performs over-representation analysis.
###############################################################################

# 01. Packages and data #######################################################

library(tidyverse)
library(seqwrap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(marginaleffects)


library(cowplot)
library(ggtext)
library(ComplexUpset)
source("figures/figure-opts.R")

models_data <- readRDS("data-out/pillon-models.RDS")


# Count number of successful models
summaries <- lapply(models_data[1:5], seqwrap_summarise, verbose = FALSE)

# Extract successful models
targets_models <- bind_rows(
  summaries$m1_results$evaluations |>
    mutate(model = "m1"),
  summaries$m2_results$evaluations |>
    mutate(model = "m2"),
  summaries$m3_results$evaluations |>
    mutate(model = "m3"),
  summaries$m4_results$evaluations |>
    mutate(model = "m4"),
  summaries$m5_results$evaluations |>
    mutate(model = "m5")
) |>

  # Filter singular fits and bad convergence
  filter(!(model == "m3" & isSingular == TRUE)) |>
  # Poisson and NB models do not have pdHess == FALSE or convergence != 0
  # but for failsafe
  filter(
    !(model %in%
      c("m1", "m2", "m4", "m5") &
      (pdHess == FALSE | convergence != 0))
  ) |>
  group_by(model) |>
  distinct(target)

# Sums number of successful estimated targets  #
n_estimates <- targets_models |>
  ungroup() |>
  summarise(.by = model, n = n()) |>
  # Add total
  mutate(
    total = nrow(models_data$filtered_counts),
    estimated_prop = 100 * (n / total),
    non_est_prop = (100 - estimated_prop)
  )


saveRDS(n_estimates, "data-out/non-est-prop.RDS")

# 01. Plot A - B number of estimates and DEGs ###############################

p1 <- n_estimates |>
  mutate(
    non_est = paste0("n = ", total - n),
    model = factor(
      model,
      levels = c("m1", "m2", "m3", "m4", "m5"),
      labels = c(
        "Negative binomial\n(non-informed)",
        "Regularized\nNegative binomial",
        "Gaussian transformed\ncounts",
        "Poisson OLRE\n(non-informed)",
        "Regularized\nPoisson OLRE"
      )
    )
  ) |>
  ggplot(aes(model, non_est_prop, fill = model)) +
  geom_bar(stat = "identity", width = 0.4) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +

  geom_text(aes(model, non_est_prop + 1.1, label = non_est), size = 2.5) +

  theme_classic(base_size = 8) +

  scale_fill_manual(
    values = c(colors[1], colors[2], colors[4], colors[1], colors[2])
  ) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(y = "Missing estimates (%)") +
  coord_flip()


## AIC per model

evaluations <- bind_rows(
  summaries$m1_results$evaluations |>
    mutate(model = "m1"),
  summaries$m2_results$evaluations |>
    mutate(model = "m2"),

  summaries$m4_results$evaluations |>
    mutate(model = "m4"),
  summaries$m5_results$evaluations |>
    mutate(model = "m5")
)


# Model AIC and compare across targets/models

# Calculate how many targets have lower AIC when comparing models

delta_aic <- evaluations |>
  dplyr::select(target, model, aic) |>
  pivot_wider(names_from = model, values_from = aic) |>

  mutate(
    m4_m1 = if_else(m4 - m1 < 0, TRUE, FALSE),
    m5_m2 = if_else(m5 - m2 < 0, TRUE, FALSE)
  ) |>
  dplyr::select(target, m4_m1, m5_m2) |>
  pivot_longer(cols = m4_m1:m5_m2) |>
  summarise(.by = name, n = sum(value, na.rm = TRUE)) |>
  mutate(
    lab = c("NB > POLRE: ", "NB-R > POLRE-R: "),
    n_total = nrow(models_data$filtered_counts),
    prop = 100 * (n / n_total),
    prop = paste0(lab, round(prop, 0), "%"),
    xmin = c(1, 2),
    xmax = c(3, 4),
    ycoord = c(1065, 1062)
  )


aic_mod <- lme4::lmer(aic ~ model + (1 | target), data = evaluations)

marginaleffects::avg_comparisons(aic_mod)

pred <- marginaleffects::avg_predictions(aic_mod, re.form = NULL, by = "model")


p1b <- data.frame(pred) |>

  mutate(
    model = factor(
      model,
      levels = c("m1", "m2", "m4", "m5"),
      labels = c(
        "Negative binomial\n(NB)",
        "Regularized Negative\nbinomial (NB-R)",

        "Poisson OLRE\n(POLRE)",
        "Regularized Poisson\nOLRE (POLRE-R)"
      )
    )
  ) |>

  ggplot(aes(estimate, model)) +

  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.2
  ) +

  geom_point(aes(shape = model, fill = model), size = 3) +

  theme_classic(base_size = 8) +

  scale_fill_manual(values = c(colors[1], colors[2], colors[1], colors[2])) +

  scale_shape_manual(values = c(22, 22, 23, 25, 25)) +

  scale_x_continuous(limits = c(1040, 1080), expand = c(0, 0)) +

  geom_segment(
    data = delta_aic,
    aes(y = xmin, yend = xmax, x = ycoord, xend = ycoord)
  ) +
  geom_segment(
    data = delta_aic,
    aes(y = xmin, yend = xmin, x = ycoord, xend = ycoord - 1)
  ) +
  geom_segment(
    data = delta_aic,
    aes(y = xmax, yend = xmax, x = ycoord, xend = ycoord - 1)
  ) +

  geom_text(
    data = delta_aic,
    aes(
      y = xmax - (xmax - xmin) / 2 + 0.5,
      x = ycoord + 0.5,
      hjust = 0,
      label = prop
    ),
    size = 2.5
  ) +

  labs(x = "AIC (Mean &pm; CI) ") +

  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(),
    legend.position = "none"
  )


# Summaries for gene set analysis and sig genes ################################

stat <- bind_rows(
  summaries$m1_results$summaries |>
    mutate(model = "m1"),
  summaries$m2_results$summaries |>
    mutate(model = "m2"),
  summaries$m3_results$summaries |>
    mutate(model = "m3"),
  summaries$m4_results$summaries |>
    mutate(model = "m4"),
  summaries$m5_results$summaries |>
    mutate(model = "m5")
) |>

  # This filter away singular fits from m3
  inner_join(targets_models) |>

  filter(term %in% c("groupT2D", "timepost:groupT2D", "timerec:groupT2D")) |>

  # Adjust p-values and calculate MSD
  mutate(.by = c(model, term), fdr = p.adjust(p.value, method = "fdr")) |>
  # Approximate CI
  mutate(
    lwr = estimate - qnorm(0.975) * std.error,
    upr = estimate + qnorm(0.972) * std.error,
    msd = if_else(estimate > 0, lwr, -upr)
  )


saveRDS(stat, "data-out/estimates-pillon.RDS")


# 00. Number of significant genes per model/term upset plot

sig_df <- stat |>
  filter(fdr < 0.05) |>
  mutate(
    model = factor(
      model,
      levels = c("m1", "m2", "m3", "m4", "m5"),
      labels = c(
        "Negative binomial (non-informed)",
        "Regularized Negative binomial",
        "Gaussian transformed counts",
        "Poisson OLRE (non-informed)",
        "Regularized Poisson OLRE"
      )
    )
  )


sig_genes_df <- stat |>
  filter(fdr < 0.05) |>
  mutate(
    model = factor(
      model,
      levels = c("m1", "m2", "m3", "m4", "m5"),
      labels = c(
        "Negative binomial (non-informed)",
        "Regularized Negative binomial",
        "Gaussian transformed counts",
        "Poisson OLRE (non-informed)",
        "Regularized Poisson OLRE"
      )
    )
  ) |>
  dplyr::select(target, model) |>
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = model,
    values_from = value,
    values_fill = FALSE
  )


p2 <- upset(
  sig_genes_df,
  intersect = setdiff(names(sig_genes_df), "target"),
  name = "",
  height_ratio = 1.2,
  width_ratio = 0.2,
  stripes = c(
    alpha(colors[4], 0.2),
    alpha(colors[1], 0.2),
    alpha(colors[1], 0.2),
    alpha(colors[2], 0.2),
    alpha(colors[2], 0.2)
  ),

  set_sizes = (upset_set_size(
    geom = geom_bar(width = 0.4),
    position = "right"
  ) +
    geom_text(
      aes(label = after_stat(count)),
      hjust = -0.2,
      stat = 'count',
      size = 2.5
    ) +
    expand_limits(y = 4000) +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())),

  matrix = (intersection_matrix() +
    theme(axis.text.y = element_text(color = "black"))),

  base_annotations = list(
    'Intersection size' = (intersection_size(
      text_mapping = aes(
        label = !!upset_text_percentage()
      ),
      text = list(size = 2.5)
    ) +
      ylim(c(0, 1200)) +
      ylab('') +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 8, color = "black")
      ))
  ),

  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(text = element_text(size = 8))
    )
  ),
  min_size = 10
)


# 03. Performs over-representation analysis. ################################
# ora <- function(Term = "timerec:groupT2D", Model = "m5") {
#
#   # Pull genes with fdr < 0.01
#   gene_list <- stat |>
#     filter(term == Term, model == Model) |>
#     filter(fdr < 0.05) |>
#     pull(target)
#
#   # Pull universe
#   universe <- stat |>
#     filter(term == Term, model == Model) |>
#     pull(target)
#
#   gene_list <- bitr(gene_list,
#                     fromType = "SYMBOL",
#                     toType = "ENSEMBL",
#                     OrgDb = org.Hs.eg.db)
#   universe <- bitr(universe,
#                    fromType = "SYMBOL",
#                    toType = "ENSEMBL",
#                    OrgDb = org.Hs.eg.db)
#
#
#
#   if (length(gene_list) > 0) {
#
#     ora_results <-   enrichGO(gene         = gene_list$ENSEMBL,
#                               universe     = universe$ENSEMBL,
#                               OrgDb         = org.Hs.eg.db,
#                               keyType       = 'ENSEMBL',
#                               ont           = "BP",
#                               pAdjustMethod = "BH",
#                               pvalueCutoff  = 0.01,
#                               qvalueCutoff  = 0.05)
#
#   } else { ora_results <- NULL }
#
#
#
# return(ora_results)
#
#
# }
#
#
#
#
#
# term_model <- expand_grid(term = "timerec:groupT2D",
#                           model = unique(stat$model))
#
# results <- list()
# for(i in seq_along(1:nrow(term_model))) {
#
#
#   res_temp <- ora(Term = term_model[i, 1][[1]],
#                   Model = term_model[i, 2][[1]])
#
#   if (!is.null(res_temp)) {
#     results[[i]] <- res_temp@result |>
#       data.frame(row.names = NULL) |>
#       mutate(term = term_model[i, 1][[1]],
#              model = term_model[i, 2][[1]])
#   } else {
#
#     results[[i]] <- data.frame(term = term_model[i, 1][[1]],
#                                 model = term_model[i, 2][[1]])
#
#   }
#
#
#   print(paste0("Iter ", i, " of ", nrow(term_model)))
#
#
# }
#
#
# # Number of significant gene sets in ORA per model
# ora_terms_df <- bind_rows(results) |>
#   filter(p.adjust < 0.01) |>
#   dplyr::select(model, ID) |>
#   print()
#
# ora_terms_list <- split(ora_terms_df$ID,
#       ora_terms_df$model)
#
#
# upset(fromList(ora_terms_list), order.by = "freq")
#
#
# bind_rows(results) |>
#   # Extract common important terms
#   mutate(.by = c(ID),
#          effect_size = mean(-log10(p.adjust))) |>
#   mutate(.by = model,
#          rank = rank(-(-log10(pvalue)))) |>
#
#   dplyr::select(ID:Description, zScore:qvalue, Count, term, model, effect_size, rank) |>
#   group_by(model) |>
#   slice_max(order_by = -rank, n = 4) |>
#
#   ggplot(aes(zScore, ID, color = model)) + geom_point()
#
#

# 04. GSEA Analaysis #######################################################

gsea <- function(Term = "timerec:groupT2D", Model = "m2") {
  gene_list <- stat |>
    filter(term == Term, model == Model) |>
    arrange(-estimate) |>
    dplyr::select(target, estimate)

  geneList <- gene_list$estimate
  names(geneList) <- as.character(gene_list$target)

  gse_res <- gseGO(
    geneList = geneList,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )

  df <- gse_res@result |>
    mutate(model = Model, term = Term)

  return(list(results = gse_res, df = df))
}


res1 <- gsea(Term = "timerec:groupT2D", Model = "m1")
res2 <- gsea(Term = "timerec:groupT2D", Model = "m2")
res3 <- gsea(Term = "timerec:groupT2D", Model = "m3")
res4 <- gsea(Term = "timerec:groupT2D", Model = "m4")
res5 <- gsea(Term = "timerec:groupT2D", Model = "m5")

# Upset plot of gsea with fdr < alpha

gsea_results <- bind_rows(
  list(
    data.frame(res1$df, row.names = NULL),
    data.frame(res2$df, row.names = NULL),
    data.frame(res3$df, row.names = NULL),
    data.frame(res4$df, row.names = NULL),
    data.frame(res5$df, row.names = NULL)
  )
)


saveRDS(gsea_results, "data-out/gsea_results.RDS")


gsea_df <- bind_rows(
  list(
    data.frame(res1$df, row.names = NULL),
    data.frame(res2$df, row.names = NULL),
    data.frame(res3$df, row.names = NULL),
    data.frame(res4$df, row.names = NULL),
    data.frame(res5$df, row.names = NULL)
  )
) |>
  filter(p.adjust < 0.01) |>
  dplyr::select(model, ID) |>
  print()


## Create axis with fold changes (average over all models)

anno_df <- stat |>
  filter(term == "timerec:groupT2D") |>
  summarise(.by = target, fc = mean(estimate)) |>
  mutate(
    rank = rank(-fc),
    distance.neg = abs(fc - (-0.5)),
    distance.mid = abs(fc - 0),
    distance.pos = abs(fc - 0.5),
    axis.neg = if_else(distance.neg == min(distance.neg), rank, NA),
    axis.mid = if_else(distance.mid == min(distance.mid), rank, NA),
    axis.pos = if_else(distance.pos == min(distance.pos), rank, NA)
  ) |>
  filter(
    distance.neg == min(distance.neg) |
      distance.pos == min(distance.pos) |
      distance.mid == min(distance.mid)
  )


fc_axis <- stat |>
  filter(term == "timerec:groupT2D") |>
  summarise(.by = target, fc = mean(estimate)) |>
  mutate(rank = rank(-fc)) |>

  ggplot(aes(rank, y = 1)) +
  geom_segment(aes(x = rank, xend = rank, y = 1.05, yend = 1.95, color = fc)) +

  annotate(
    "segment",
    x = anno_df$axis.neg[[3]],
    xend = anno_df$axis.pos[[1]],
    y = 1,
    yend = 1
  ) +
  # high
  annotate(
    "segment",
    x = anno_df$axis.pos[[1]],
    xend = anno_df$axis.pos[[1]],
    y = 1,
    yend = 0.8
  ) +
  # mid
  annotate(
    "segment",
    x = anno_df$axis.mid[[2]],
    xend = anno_df$axis.mid[[2]],
    y = 1,
    yend = 0.8
  ) +
  annotate(
    "segment",
    x = anno_df$axis.neg[[3]],
    xend = anno_df$axis.neg[[3]],
    y = 1,
    yend = 0.8
  ) +

  annotate(
    "text",
    x = c(
      anno_df$axis.pos[[1]],
      anno_df$axis.mid[[2]],
      anno_df$axis.neg[[3]],
      anno_df$axis.mid[[2]]
    ),
    y = c(rep(0.5, 3), 0.15),
    size = 2.5,
    label = c("0.5", "0", "-0.5", "Fold-change")
  ) +

  scale_color_gradient2() +
  scale_y_continuous(limits = c(0, 2)) +

  theme_void() +
  theme(legend.position = "none")


# 0030199 Collagen Fibril Organization
pdat1_up <- gseaplot(res1$results, geneSetID = "GO:0030199")
pdat2_up <- gseaplot(res2$results, geneSetID = "GO:0030199")
pdat3_up <- gseaplot(res3$results, geneSetID = "GO:0030199")
pdat4_up <- gseaplot(res4$results, geneSetID = "GO:0030199")
pdat5_up <- gseaplot(res5$results, geneSetID = "GO:0030199")

##  Aerobic Electron Transport Chain
pdat1 <- gseaplot(res1$results, geneSetID = "GO:0019646")
pdat2 <- gseaplot(res2$results, geneSetID = "GO:0019646")
pdat3 <- gseaplot(res3$results, geneSetID = "GO:0019646")
pdat4 <- gseaplot(res4$results, geneSetID = "GO:0019646")
pdat5 <- gseaplot(res5$results, geneSetID = "GO:0019646")


gsea_plot_up <- bind_rows(
  pdat1_up[[1]]$data |>
    mutate(model = "m1"),
  pdat2_up[[1]]$data |>
    mutate(model = "m2"),
  pdat3_up[[1]]$data |>
    mutate(model = "m3"),
  pdat4_up[[1]]$data |>
    mutate(model = "m4"),
  pdat5_up[[1]]$data |>
    mutate(model = "m5")
) |>
  mutate(
    model = factor(
      model,
      levels = c("m1", "m2", "m3", "m4", "m5"),
      labels = c(
        "Negative binomial (non-informed)",
        "Regularized Negative binomial",
        "Gaussian transformed counts",
        "Poisson OLRE (non-informed)",
        "Regularized Poisson OLRE"
      )
    )
  ) |>

  mutate(.by = model, xrelative = x / max(x)) |>

  ggplot(aes(xrelative, runningScore, color = model)) +
  geom_line() +
  theme_classic(base_size = 8) +

  scale_color_manual(values = c(colors)) +

  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0, r = 2),
    legend.margin = margin(r = 20),
    legend.title = element_blank()
  ) +
  labs(
    y = "GSEA running score",
    subtitle = "Gene ontology: Collagen Fibril Organization",
    color = "Model"
  )


gsea_plot_down <- bind_rows(
  pdat1[[1]]$data |>
    mutate(model = "m1"),
  pdat2[[1]]$data |>
    mutate(model = "m2"),
  pdat3[[1]]$data |>
    mutate(model = "m3"),
  pdat4[[1]]$data |>
    mutate(model = "m4"),
  pdat5[[1]]$data |>
    mutate(model = "m5")
) |>
  mutate(
    model = factor(
      model,
      levels = c("m1", "m2", "m3", "m4", "m5"),
      labels = c(
        "Negative binomial (non-informed)",
        "Regularized Negative binomial",
        "Gaussian transformed counts",
        "Poisson OLRE (non-informed)",
        "Regularized Poisson OLRE"
      )
    )
  ) |>

  mutate(.by = model, xrelative = x / max(x)) |>

  ggplot(aes(xrelative, runningScore, color = model)) +
  geom_line() +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0),
    legend.position = "none"
  ) +

  scale_color_manual(values = c(colors)) +

  labs(
    y = "GSEA running score",
    subtitle = "Gene ontology: Aerobic Electron Transport Chain",
    color = "Model"
  )


p3 <- plot_grid(
  NULL,
  plot_grid(
    gsea_plot_up,
    fc_axis,
    NULL,
    gsea_plot_down,
    ncol = 1,
    rel_heights = c(1, 0.6, 0.05, 1),
    align = "vh"
  ),
  NULL,
  rel_widths = c(0.05, 1, 0.05),
  ncol = 3
)


## Combine plot #############################

figure4 <- plot_grid(
  plot_grid(NULL, p1, p1b, NULL, ncol = 4, rel_widths = c(0.05, 0.7, 1, 0.05)),
  NULL,
  p2,
  p3,
  NULL,
  rel_heights = c(0.6, 0.1, 1, 1, 0.05),

  ncol = 1
) +
  annotate(
    "text",
    x = c(0.03, 0.4, 0.05, 0.03),
    y = c(0.98, 0.98, 0.7, 0.42),
    label = c("A", "B", "C", "D")
  )


## Save RDS ##
saveRDS(figure4, "figures/figure-4.RDS", compress = "xz")


## Save PDF ##
ggsave(
  "figures/figure-4.pdf",
  figure4,
  device = cairo_pdf,
  height = 200,
  width = 170,
  units = "mm"
)

# MA-PLOT ####################################################

#
#
#
#  estimates <- bind_rows(
#    summaries$m1_results$summaries |>
#      mutate(model = "m1"),
#    summaries$m2_results$summaries |>
#      mutate(model = "m2"),
#
#    summaries$m4_results$summaries |>
#      mutate(model = "m4"),
#    summaries$m5_results$summaries |>
#      mutate(model = "m5") )
#
#
#  estimates |>
#    filter(term != "sd__(Intercept)") |>
#
#     dplyr::select(target, model, term, estimate) |>
#    pivot_wider(names_from = term, values_from = estimate) |>
#
#    # Calculate mean at recovery from exercise #
#    mutate(mean.timerec.ngt = exp(`(Intercept)` + timerec),
#           mean.timerec.t2d = exp(`(Intercept)` + timerec + groupT2D + `timerec:groupT2D`),
#           mean.timerec = log((mean.timerec.ngt + mean.timerec.t2d) / 2)) |>
#    dplyr::select(target, model, mean.timerec) |>
#    inner_join(
#      estimates |>
#                   filter(term == "timerec:groupT2D") |>
#                   mutate(fdr = p.adjust(p.value, method = "fdr")) |>
#
#                   dplyr::select(target, model, estimate, fdr)
#      )  |>
#    mutate(sig = if_else(fdr < 0.05, "s", "ns") ) |>
#
#    filter(model %in% c("m2", "m5")) |>
#    dplyr::select(-mean.timerec) |>
#
#    pivot_wider(names_from = model,
#                values_from = c(estimate, fdr, sig)) |>
#    mutate(concl = if_else(sig_m2 == sig_m5, "same", "diff")) |>
#
#
#
#
#    ggplot(aes(-log10(fdr_m2), -log10(fdr_m5), color = concl)) + geom_point(alpha = 0.2)
#
#
#    ggplot(aes(mean.timerec, estimate, color = sig)) +geom_point() + facet_grid(model ~ .)
#
#
#
# temp <-    estimates |>
#      filter(term != "sd__(Intercept)") |>
#
#      dplyr::select(target, model, term, estimate) |>
#      pivot_wider(names_from = term, values_from = estimate) |>
#
#      # Calculate mean at recovery from exercise #
#      mutate(mean.timerec.ngt = exp(`(Intercept)` + timerec),
#             mean.timerec.t2d = exp(`(Intercept)` + timerec + groupT2D + `timerec:groupT2D`),
#             mean.timerec = log((mean.timerec.ngt + mean.timerec.t2d) / 2)) |>
#      dplyr::select(target, model, mean.timerec) |>
#      inner_join(
#        estimates |>
#          filter(term == "timerec:groupT2D") |>
#          mutate(fdr = p.adjust(p.value, method = "fdr")) |>
#
#          dplyr::select(target, model, estimate, fdr)
#      )  |>
#      mutate(sig = if_else(fdr < 0.05, "s", "ns") ) |>
#
#      filter(model %in% c("m2", "m5")) |>
#
#      pivot_wider(names_from = model,
#                  values_from = c(mean.timerec, estimate, fdr, sig)) |>
#      mutate(concl = if_else(sig_m2 == sig_m5, "same", "diff")) |>
#      print()
#
#
#
#
#
#
# temp |>
#
#      ggplot(aes(mean.timerec_m2, estimate_m2)) +
#   geom_point(alpha = 0.2) +
#   geom_point(data = filter(temp, concl == "diff"),
#              aes(mean.timerec_m5, estimate_m5), size = 3, color = "red") +
# geom_point(data = filter(temp, concl == "diff"),
#            aes(mean.timerec_m2, estimate_m2), size = 3, color = "blue") +
#   geom_segment(data = filter(temp, concl == "diff"),
#              aes(x = mean.timerec_m2, xend = mean.timerec_m5,
#                  y = estimate_m2, yend = estimate_m5))
#
#
#
#
#
#
# bind_rows(
#   summaries$m1_results$evaluations |>
#     mutate(model = "m1"),
#   summaries$m2_results$evaluations |>
#     mutate(model = "m2"),
#
#   summaries$m4_results$evaluations |>
#     mutate(model = "m4"),
#   summaries$m5_results$evaluations |>
#     mutate(model = "m5") )  |>
#
#   filter(model %in% c("m2", "m5")) |>
#   dplyr::select(model, target:dispersion.se) |>
#   pivot_wider(names_from = model,
#               values_from = c(dispersion, dispersion.se)) |>
#
#   mutate(disp.diff = dispersion_m2 - dispersion_m5) |>
#
#   inner_join(temp) |>
#
#
#
#   ggplot(aes(dispersion_m2, dispersion_m5, color = concl)) + geom_point(alpha = 0.2)
#
#   print()
#
#
#
#obs <-  data.frame(target = models_data$filtered_counts[,1],
#            mean = rowMeans(models_data$filtered_counts[,-1]),
#            var.obs = apply(models_data$filtered_counts[,-1], 1, var))
#
#
## Extract SD of ran effects
#
#pois_sdid <- estimates |>
#  dplyr::filter(model %in% c("m4", "m5"),
#                      term == "sd__(Intercept)",
#                      group == "id") |>
#  dplyr::select(target, model, sd.id = estimate) |>
#  tidyr::complete(model = c("m1", "m2"),
#                  target = unique(estimates$target)) |>
#  print()
#
#
#
#bind_rows(
#  summaries$m1_results$evaluations |>
#    mutate(model = "m1"),
#  summaries$m2_results$evaluations |>
#    mutate(model = "m2"),
#
#  summaries$m4_results$evaluations |>
#    mutate(model = "m4"),
#  summaries$m5_results$evaluations |>
#    mutate(model = "m5") ) |>
#
#  inner_join(pois_sdid) |>
#
#  mutate(var = if_else(is.na(olre.sd),
#                       exp(log_mu) + exp(log_mu)^2 / exp(dispersion),
#                       exp(log_mu) + exp(log_mu)^2 * (olre.sd^2 + sd.id^2))) |>
#  dplyr::select(target, model, log_mu, var) |>
#  inner_join(obs) |>
#
##  summarise(.by = model,
##            c = cor(log(var), log(var.obs)))
#
#
#  ggplot(aes(log(var.obs), log(var))) + geom_point() +
#  facet_wrap(~ model) +
#
#  geom_abline(slope = 1, intercept = 1)
#
#
#  print()
#
#
#
#
#
#
#  bind_rows(
#    summaries$m1_results$evaluations |>
#      mutate(model = "m1"),
#    summaries$m2_results$evaluations |>
#      mutate(model = "m2"),
#
#    summaries$m4_results$evaluations |>
#      mutate(model = "m4"),
#    summaries$m5_results$evaluations |>
#      mutate(model = "m5") ) |>
#
#    inner_join(pois_sdid) |>
#
#    mutate(var = if_else(is.na(olre.sd),
#                         exp(log_mu) + exp(log_mu)^2 / exp(dispersion),
#                         exp(log_mu) + exp(log_mu)^2 * (olre.sd^2 + sd.id^2))) |>
#
#    mutate(dispersion = if_else(model %in% c("m1", "m2"),
#                                dispersion,
#                                olre.sd^2 + sd.id^2) ) |>
#
#    dplyr::select(target, model, dispersion) |>
#
#    pivot_wider(names_from = model, values_from = dispersion) |>
#
#
# #   summarise(c = cor(-m2, log(m5), use = "complete.obs"))
#
#    ggplot(aes(-m2, log(m5))) + geom_point()
#
#    dplyr::select()
#
#
#
#
#
#
#estimates |>
#  filter(model %in% c("m5", "m4")) |>
#  filter(group == "seq_sample_id") |>
#  dplyr::select(model, target, estimate) |>
#  inner_join(obs) |>
#  mutate(polre = mean + estimate^2) |>
#
#  ggplot(aes(log(mean), log(polre))) + geom_point() +
#  facet_wrap(~model)
#
#
#  summarise(.by = model,
#            c = cor(log(var), log(polre)))
#
#  ggplot(aes(log(var), log(polre))) + geom_point() +
#  facet_wrap(~model)
#
#
#bind_rows(
#  summaries$m1_results$evaluations |>
#    mutate(model = "m1"),
#  summaries$m2_results$evaluations |>
#    mutate(model = "m2")) |>
#  dplyr::select(model, target, dispersion) |>
#  inner_join(obs) |>
#  mutate(nb_var = mean + mean^2 / exp(dispersion)) |>
#
#  ggplot(aes(log(mean), log(nb_var))) + geom_point() +
#  facet_wrap(~model)
#
#
#  summarise(.by = model,
#            c = cor(log(var), log(nb_var)))
#
#
#  ggplot(aes(log(var), log(nb_var))) + geom_point() +
#  facet_wrap(~model)
#
