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
library(ggupset)
source(here::here("analysis/figures/figure-opts.R"))


# Check that all models exists
models <- paste0("m", 1:5, "_results.RDS")

if (!all(models %in% list.files(here::here("analysis/data/derived_data/")))) {
  stop(
    "Not all models have results stored in 'analysis/data/derived_data/'.
           Re-run models using 'analysis/paper/make-docs.R'"
  )
}


# Collect all model results
models_data <- list()
for (i in 1:length(models)) {
  models_data[[i]] <- readRDS(
    here::here(
      paste0("analysis/data/derived_data/", models[i])
    )
  )
  names(models_data)[i] <- paste0("m", i, "_results")
}


# Get number of filtered counts
filtered_counts <- readRDS(here::here(
  "analysis/data/derived_data/filtered_counts.RDS"
))


# Count number of successful models
summaries <- lapply(models_data[1:5], seqwrap_summarise, verbose = FALSE)


# Extract successful models
targets_models <- bind_rows(
  summaries$m1_results$evaluations |>
    mutate(model = "m1"),
  summaries$m2_results$evaluations |>
    mutate(model = "m2"),

  summaries$m4_results$evaluations |>
    mutate(model = "m4"),
  summaries$m5_results$evaluations |>
    mutate(model = "m5")
) |>
  filter(convergence == 0 & pdHess == TRUE) |>
  group_by(model) |>
  distinct(target) |>
  ungroup() |>

  # For model m3, is singular is the main culprit
  bind_rows(
    summaries$m3_results$evaluations |>
      mutate(model = "m3") |>

      filter(isSingular == "FALSE") |>
      group_by(model) |>
      distinct(target) |>
      ungroup()
  )


# Sums number of successful estimated targets  #
n_estimates <- targets_models |>
  ungroup() |>
  summarise(.by = model, n = n()) |>
  # Add total
  mutate(
    total = nrow(filtered_counts),
    estimated_prop = 100 * (n / total),
    non_est_prop = (100 - estimated_prop)
  )


saveRDS(n_estimates, here::here("analysis/data/derived_data/n_estimates.RDS"))

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
    n_total = nrow(filtered_counts),
    prop = 100 * (n / n_total),
    prop = paste0(lab, round(prop, 0), "%"),
    xmin = c(1, 2),
    xmax = c(3, 4),
    ycoord = c(1065, 1062)
  )


aic_mod <- lme4::lmer(aic ~ model + (1 | target), data = evaluations)


pred <- marginaleffects::avg_predictions(aic_mod, re.form = NA, by = "model")


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


saveRDS(stat, here::here("analysis/data/derived_data/estimates-pillon.RDS"))


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


# Calculate sum of all declared significant genes
total_n <- sig_genes_df |>
  pivot_longer(cols = -target) |>
  filter(value) |>
  distinct(target) |>
  nrow()


upset <- sig_genes_df |>
  pivot_longer(cols = -target) |>
  filter(value) |>
  group_by(target) |>
  summarise(
    model = list(name)
  ) |>

  ggplot(aes(x = model)) +
  geom_bar() +
  scale_x_upset() +
  geom_text(
    stat = 'count',
    aes(
      label = paste0(
        round(100 * (after_stat(count) / total_n), 0),
        "%"
      )
    ),
    vjust = -1,
    size = 2.5
  ) +
  theme(
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 1200),
    expand = c(0, 0),
    breaks = c(0, 250, 500, 750, 1000)
  )

## Set size figure

set_size <- sig_genes_df |>
  pivot_longer(cols = -target) |>
  filter(value) |>
  summarise(.by = name, setsize = sum(value)) |>

  ggplot(aes(x = setsize, y = name)) +
  geom_bar(stat = "identity", width = 0.6, aes(fill = name)) +
  geom_text(
    position = position_nudge(x = 450),
    aes(label = setsize),
    size = 3
  ) +
  scale_x_continuous(limits = c(0, 4000), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      alpha(colors[4], 0.7),
      alpha(colors[1], 0.7),
      alpha(colors[1], 0.7),
      alpha(colors[2], 0.7),
      alpha(colors[2], 0.7)
    )
  ) +
  theme(
    axis.title = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none"
  )


p2 <- plot_grid(
  upset,
  plot_grid(NULL, set_size, NULL, rel_heights = c(1, 0.62, 0.01), ncol = 1),
  rel_widths = c(1, 0.3),
  ncol = 2
)


# 04. GSEA Analaysis #######################################################

gsea <- function(Term = "timerec:groupT2D", Model = "m5") {
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


if (!file.exists(here::here("analysis/data/derived_data/gsea_results.RDS"))) {
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

  gsea_list <- list(
    res1 = res1,
    res2 = res2,
    res3 = res3,
    res4 = res4,
    res5 = res5
  )

  saveRDS(gsea_list, here::here("analysis/data/derived_data/gsea_list.RDS"))
  saveRDS(
    gsea_results,
    here::here("analysis/data/derived_data/gsea_results.RDS")
  )
}

gsea_list <- readRDS(here::here("analysis/data/derived_data/gsea_list.RDS"))
gsea_results <- readRDS(here::here(
  "analysis/data/derived_data/gsea_results.RDS"
))


gsea_df <- gsea_results |>
  filter(p.adjust < 0.01) |>
  dplyr::select(model, ID)


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
pdat1_up <- gseaplot(gsea_list$res1$results, geneSetID = "GO:0030199")
pdat2_up <- gseaplot(gsea_list$res2$results, geneSetID = "GO:0030199")
pdat3_up <- gseaplot(gsea_list$res3$results, geneSetID = "GO:0030199")
pdat4_up <- gseaplot(gsea_list$res4$results, geneSetID = "GO:0030199")
pdat5_up <- gseaplot(gsea_list$res5$results, geneSetID = "GO:0030199")

##  Aerobic Electron Transport Chain
pdat1 <- gseaplot(gsea_list$res1$results, geneSetID = "GO:0019646")
pdat2 <- gseaplot(gsea_list$res2$results, geneSetID = "GO:0019646")
pdat3 <- gseaplot(gsea_list$res3$results, geneSetID = "GO:0019646")
pdat4 <- gseaplot(gsea_list$res4$results, geneSetID = "GO:0019646")
pdat5 <- gseaplot(gsea_list$res5$results, geneSetID = "GO:0019646")


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
    y = c(0.98, 0.98, 0.7, 0.38),
    label = c("A", "B", "C", "D")
  )


## Save RDS ##
saveRDS(figure4, here::here("analysis/figures/figure-4.RDS"), compress = "xz")


## Save PDF ##
ggsave(
  here::here("analysis/figures/figure-4.pdf"),
  figure4,
  device = cairo_pdf,
  height = 200,
  width = 170,
  units = "mm"
)

ggsave(
      here::here("analysis/figures/figure-4.png"),
      figure4,
      bg = "white",
      dpi = 600,
      height = 200,
      width = 170,
      units = "mm"
)

