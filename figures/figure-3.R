
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
source("figures/figure-opts.R")

models_data  <-  readRDS("data-out/pillon-models.RDS")


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
    mutate(model = "m5") ) |>

  # Filter singular fits and bad convergence
  filter(!(model == "m3" & isSingular == TRUE)) |>
  # Poisson and NB models do not have pdHess == FALSE or convergence != 0
  # but for failsafe
  filter(!(model %in% c("m1", "m2", "m4", "m5") &
             (pdHess == FALSE | convergence != 0))) |>
  group_by(model) |>
  distinct(target)

# Sums number of successful estimated targets  #
n_estimartes <- targets_models |>
  ungroup() |>
  summarise(.by = model,
            n = n()) |>
  # Add total
  mutate(total = nrow(models_data$filtered_counts),
         estimated_prop = 100 * (n / total),
         non_est_prop = (100 - estimated_prop))



p1 <- n_estimartes |>
  mutate(non_est = paste0("n = ", total - n),
  model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                 labels = c("Negative binomial (non-informed)",
                            "Regularized Negative binomial",
                            "Gaussian transformed counts",
                            "Poisson OLRE (non-informed)",
                            "Regularized Poisson OLRE"))) |>
  ggplot(aes(model, non_est_prop, fill = model)) +
  geom_bar(stat = "identity", width = 0.4) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +

  geom_text(aes(model, non_est_prop + 0.5, label = non_est),
            size = 3) +

  theme_classic() +

  scale_fill_manual(values = c(colors[1], colors[2], colors[4], colors[1], colors[2])) +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
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
    mutate(model = "m5") )


# Model AIC and compare across targets/models

# Calculate how many targets have lower AIC when comparing models

delta_aic <- evaluations |>
  dplyr::select(target, model, aic) |>
  pivot_wider(names_from = model,
              values_from = aic) |>

  mutate(
         m4_m1 = if_else(m4 - m1 < 0, TRUE, FALSE),
         m5_m2 = if_else(m5 - m2 < 0, TRUE, FALSE)) |>
  dplyr::select(target, m4_m1, m5_m2) |>
  pivot_longer(cols = m4_m1:m5_m2) |>
  summarise(.by = name,
            n = sum(value, na.rm = TRUE)) |>
  mutate(lab = c("Model 4 < Model 1: ",
                 "Model 5 < Model 2: "),
               n_total = nrow(models_data$filtered_counts),
         prop = 100 * (n / n_total),
         prop = paste0(lab, round(prop,1), "%"),
         xmin = c(1, 2 ),
         xmax = c( 3, 4),
         ycoord = c(1066, 1062))



aic_mod <- lme4::lmer(aic ~ model + (1|target),
           data = evaluations)

marginaleffects::avg_comparisons(aic_mod)

pred <- marginaleffects::avg_predictions(aic_mod, re.form = NULL,
                                 by = "model")



p1b <- data.frame(pred) |>

  mutate(model = factor(model, levels = c("m1", "m2", "m4", "m5"),
                 labels = c("Negative binomial (non-informed)",
                            "Regularized Negative binomial",

                            "Poisson OLRE (non-informed)",
                            "Regularized Poisson OLRE"))) |>


  ggplot(aes(model, estimate)) +

  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.2) +

  geom_point(aes(shape = model, fill = model),
             size = 3) +

  theme_classic() +

  scale_fill_manual(values = c(colors[1], colors[2], colors[1], colors[2])) +
  scale_color_manual(values = c("gray40", "gray60", "gray80")) +
  scale_shape_manual(values = c(22,22, 23, 25, 25)) +

  scale_y_continuous(limits = c(1040, 1070),
                     expand = c(0,0)) +

  geom_segment(data = delta_aic,
               aes(x = xmin, xend = xmax,
                   y = ycoord, yend = ycoord)) +
  geom_segment(data = delta_aic,
               aes(x = xmin, xend = xmin,
                   y = ycoord, yend = ycoord-1)) +
  geom_segment(data = delta_aic,
               aes(x = xmax, xend = xmax,
                   y = ycoord, yend = ycoord-1)) +



  geom_text(data = delta_aic,
               aes(x = xmax - (xmax - xmin)/2,
                   y = ycoord + 1,
                   label = prop)) +

  labs(y = "AIC (Mean &pm; CI) ") +

  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_markdown(),
        legend.title = element_blank())





# Summaries #
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
    mutate(model = "m5") ) |>

  # This filter away singular fits from m3
  inner_join(targets_models) |>

  filter(term %in% c("groupT2D", "timepost:groupT2D", "timerec:groupT2D")) |>

  # Adjust p-values and calculate MSD
  mutate(.by = c(model, term),
         fdr = p.adjust(p.value, method = "fdr")) |>
  # Approximate CI
  mutate(lwr = estimate - qnorm(0.975) * std.error,
         upr = estimate + qnorm(0.972) * std.error,
         msd = if_else(estimate > 0, lwr, -upr))


# 00. Number of significant genes per model/term

p2 <- stat |>


  mutate(model = factor(model, levels = c("m1", "m2", "m3", "m4", "m5"),
                        labels = c("Negative binomial (non-informed)",
                                   "Regularized Negative binomial",
                                   "Gaussian transformed counts",
                                   "Poisson OLRE (non-informed)",
                                   "Regularized Poisson OLRE"))) |>

  filter(fdr < 0.05) |>
  summarise(.by = c(model, term),
            n = n()) |>
  complete(term, model) |>
  mutate(n = if_else(is.na(n), 0, n)) |>
  filter(term == "timerec:groupT2D") |>

  arrange(model) |>
  mutate(total = nrow(models_data$filtered_counts),
         prop = 100 * (n / total),
         n_lab = paste0("n = ", n)) |>

  ggplot(aes(model, prop, fill = model)) +
  geom_bar(stat = "identity", width = 0.4) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c(colors[1], colors[2], colors[4], colors[1], colors[2])) +
  geom_text(aes(model, prop + 1, label = n_lab),
            size = 3) +

  theme_classic() +


  labs(y = "Proportion of targets with FDR < 0.05 (%)") +

  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  coord_flip()




# Comparing genes in known pathway

stat |>

  filter(target == "CCL18",
         term == "timerec:groupT2D") |>
  ggplot(aes(estimate, model)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
  geom_point()



  print()


stat |>
  filter(term == "timerec:groupT2D",
         model %in% c("m1", "m2")) |>

  dplyr::select(target, estimate, p.value, model) |>
  pivot_wider(names_from = model, values_from = c(estimate, p.value)) |>
  ggplot(aes(estimate_m1, estimate_m2)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)


stat |>
  filter(term == "timerec:groupT2D",
         model %in% c("m1", "m2")) |>

  dplyr::select(target, estimate, p.value = fdr, model) |>
  pivot_wider(names_from = model, values_from = c(estimate, p.value)) |>
  ggplot(aes(-log10(p.value_m1), -log10(p.value_m2))) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = -log10(0.05))


stat |>
  filter(term == "timerec:groupT2D",
         model %in% c("m2", "m5")) |>

  dplyr::select(target, estimate, p.value = fdr, model) |>
  pivot_wider(names_from = model, values_from = c(estimate, p.value)) |>
  ggplot(aes(-log10(p.value_m2), -log10(p.value_m5))) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = -log10(0.05))


# 03. Performs over-representation analysis. ################################


ora <- function(Term = "timerec:groupT2D", Model = "m5") {

  # Pull genes with fdr < 0.01
  gene_list <- stat |>
    filter(term == Term, model == Model) |>
    filter(fdr < 0.05) |>
    pull(target)

  # Pull universe
  universe <- stat |>
    filter(term == Term, model == Model) |>
    pull(target)

  gene_list <- bitr(gene_list,
                    fromType = "SYMBOL",
                    toType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db)
  universe <- bitr(universe,
                   fromType = "SYMBOL",
                   toType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db)



  if (length(gene_list) > 0) {

    ora_results <-   enrichGO(gene         = gene_list$ENSEMBL,
                              universe     = universe$ENSEMBL,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = 'ENSEMBL',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)

  } else { ora_results <- NULL }



return(ora_results)


}



term_model <- expand_grid(term = unique(stat$term),
                          model = unique(stat$model))

results <- list()
for(i in seq_along(1:nrow(term_model))) {


  res_temp <- ora(Term = term_model[i, 1][[1]],
                  Model = term_model[i, 2][[1]])

  if (!is.null(res_temp)) {
    results[[i]] <- res_temp@result |>
      data.frame(row.names = NULL) |>
      mutate(term = term_model[i, 1][[1]],
             model = term_model[i, 2][[1]])
  } else {

    results[[i]] <- data.frame(term = term_model[i, 1][[1]],
                                model = term_model[i, 2][[1]])

  }


  print(paste0("Iter ", i, " of ", nrow(term_model)))


}



bind_rows(results) |>

  filter(p.adjust < 0.05) |>
  summarise(.by = c(term, model),
            n_terms = n_distinct(ID)) |>

  print()



  # Extract common important terms
  mutate(.by = c(term, ID),
         effect_size = mean(-log10(p.adjust))) |>
  dplyr::select(ID:qvalue, Count, term, model, effect_size) |>

  slice_max(effect_size, n = 10) |>
  print()
  ungroup() |>

  filter(term == "timerec:groupT2D") |>
  ggplot(aes(zScore, ID, color = model)) + geom_point()

  print()

# 04. GSEA Analaysis #######################################################



gsea <- function(Term = "timerec:groupT2D", Model = "m2") {

gene_list <- stat |>
    filter(term == Term, model == Model) |>
    arrange(-estimate) |>
    dplyr::select(target, estimate)


geneList <- gene_list$estimate
names(geneList) <- as.character(gene_list$target)

gse_res <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 keyType      = "SYMBOL",
                 minGSSize    = 25,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

df <- gse_res@result |>
  mutate(model = Model,
         term = Term)

return(list(results = gse_res, df = df))

}


  res1 <- gsea(Term = "timerec:groupT2D", Model = "m1")
  res2 <- gsea(Term = "timerec:groupT2D", Model = "m2")
  res4 <- gsea(Term = "timerec:groupT2D", Model = "m4")
  res5 <- gsea(Term = "timerec:groupT2D", Model = "m5")


  ## Oxidative phosphorylation
  pdat1 <- gseaplot(res1$results, geneSetID = "GO:0006119")
  pdat2 <- gseaplot(res2$results, geneSetID = "GO:0006119")
  pdat4 <- gseaplot(res4$results, geneSetID = "GO:0006119")
  pdat5 <- gseaplot(res5$results, geneSetID = "GO:0006119")


  bind_rows(
    pdat1[[1]]$data |>
      mutate(model = "m1"),
    pdat2[[1]]$data |>
      mutate(model = "m2"),
    pdat4[[1]]$data |>
      mutate(model = "m4"),
    pdat5[[1]]$data |>
      mutate(model = "m5")

  ) |>

  ggplot(aes(x, runningScore, color = model)) + geom_line()



  ## Extracellular matrix organization
  pdat1 <- gseaplot(res1$results, geneSetID = "GO:0042110")
  pdat2 <- gseaplot(res2$results, geneSetID = "GO:0042110")
  pdat4 <- gseaplot(res4$results, geneSetID = "GO:0042110")
  pdat5 <- gseaplot(res5$results, geneSetID = "GO:0042110")


  bind_rows(
    pdat1[[1]]$data |>
      mutate(model = "m1"),
    pdat2[[1]]$data |>
      mutate(model = "m2"),
    pdat4[[1]]$data |>
      mutate(model = "m4"),
    pdat5[[1]]$data |>
      mutate(model = "m5")

  ) |>

    ggplot(aes(x, runningScore, color = model)) + geom_line()








  ## Combine plot #############################

 plot_grid(
   plot_grid(p1, p2,
            ncol = 2, rel_widths = c(1, 0.8)),
   plot_grid(p1b),
   ncol = 1)
















# MA-PLOT ####################################################














  estimates <- bind_rows(
    summaries$m1_results$summaries |>
      mutate(model = "m1"),
    summaries$m2_results$summaries |>
      mutate(model = "m2"),

    summaries$m4_results$summaries |>
      mutate(model = "m4"),
    summaries$m5_results$summaries |>
      mutate(model = "m5") )


  estimates |>
    filter(term != "sd__(Intercept)") |>

     dplyr::select(target, model, term, estimate) |>
    pivot_wider(names_from = term, values_from = estimate) |>

    # Calculate mean at recovery from exercise #
    mutate(mean.timerec.ngt = exp(`(Intercept)` + timerec),
           mean.timerec.t2d = exp(`(Intercept)` + timerec + groupT2D + `timerec:groupT2D`),
           mean.timerec = log((mean.timerec.ngt + mean.timerec.t2d) / 2)) |>
    dplyr::select(target, model, mean.timerec) |>
    inner_join(
      estimates |>
                   filter(term == "timerec:groupT2D") |>
                   mutate(fdr = p.adjust(p.value, method = "fdr")) |>

                   dplyr::select(target, model, estimate, fdr)
      )  |>
    mutate(sig = if_else(fdr < 0.05, "s", "ns") ) |>

    filter(model %in% c("m2", "m5")) |>
    dplyr::select(-mean.timerec) |>

    pivot_wider(names_from = model,
                values_from = c(estimate, fdr, sig)) |>
    mutate(concl = if_else(sig_m2 == sig_m5, "same", "diff")) |>




    ggplot(aes(-log10(fdr_m2), -log10(fdr_m5), color = concl)) + geom_point(alpha = 0.2)


    ggplot(aes(mean.timerec, estimate, color = sig)) +geom_point() + facet_grid(model ~ .)



 temp <-    estimates |>
      filter(term != "sd__(Intercept)") |>

      dplyr::select(target, model, term, estimate) |>
      pivot_wider(names_from = term, values_from = estimate) |>

      # Calculate mean at recovery from exercise #
      mutate(mean.timerec.ngt = exp(`(Intercept)` + timerec),
             mean.timerec.t2d = exp(`(Intercept)` + timerec + groupT2D + `timerec:groupT2D`),
             mean.timerec = log((mean.timerec.ngt + mean.timerec.t2d) / 2)) |>
      dplyr::select(target, model, mean.timerec) |>
      inner_join(
        estimates |>
          filter(term == "timerec:groupT2D") |>
          mutate(fdr = p.adjust(p.value, method = "fdr")) |>

          dplyr::select(target, model, estimate, fdr)
      )  |>
      mutate(sig = if_else(fdr < 0.05, "s", "ns") ) |>

      filter(model %in% c("m2", "m5")) |>

      pivot_wider(names_from = model,
                  values_from = c(mean.timerec, estimate, fdr, sig)) |>
      mutate(concl = if_else(sig_m2 == sig_m5, "same", "diff")) |>
      print()


 temp |>

      ggplot(aes(mean.timerec_m2, estimate_m2)) +
   geom_point(alpha = 0.2) +
   geom_point(data = filter(temp, concl == "diff"),
              aes(mean.timerec_m5, estimate_m5), size = 3, color = "red") +
 geom_point(data = filter(temp, concl == "diff"),
            aes(mean.timerec_m2, estimate_m2), size = 3, color = "blue") +
   geom_segment(data = filter(temp, concl == "diff"),
              aes(x = mean.timerec_m2, xend = mean.timerec_m5,
                  y = estimate_m2, yend = estimate_m5))


 bind_rows(
   summaries$m1_results$evaluations |>
     mutate(model = "m1"),
   summaries$m2_results$evaluations |>
     mutate(model = "m2"),

   summaries$m4_results$evaluations |>
     mutate(model = "m4"),
   summaries$m5_results$evaluations |>
     mutate(model = "m5") )  |>

   filter(model %in% c("m1", "m2")) |>
   dplyr::select(model, target:dispersion.se) |>
   pivot_wider(names_from = model,
               values_from = c(dispersion, dispersion.se)) |>

   mutate(disp.diff = dispersion_m1 - dispersion_m2) |>

   inner_join(temp) |>



   ggplot(aes(dispersion_m1, dispersion_m2, color = concl)) + geom_point(alpha = 0.2)

   print()



obs <-  data.frame(target = models_data$filtered_counts[,1],
            mean = rowMeans(models_data$filtered_counts[,-1]),
            var.obs = apply(models_data$filtered_counts[,-1], 1, var))


# Extract SD of ran effects

pois_sdid <- estimates |>
  dplyr::filter(model %in% c("m4", "m5"),
                      term == "sd__(Intercept)",
                      group == "id") |>
  dplyr::select(target, model, sd.id = estimate) |>
  tidyr::complete(model = c("m1", "m2"),
                  target = unique(estimates$target)) |>
  print()



bind_rows(
  summaries$m1_results$evaluations |>
    mutate(model = "m1"),
  summaries$m2_results$evaluations |>
    mutate(model = "m2"),

  summaries$m4_results$evaluations |>
    mutate(model = "m4"),
  summaries$m5_results$evaluations |>
    mutate(model = "m5") ) |>

  inner_join(pois_sdid) |>

  mutate(var = if_else(is.na(olre.sd),
                       exp(log_mu) + exp(log_mu)^2 / exp(dispersion),
                       exp(log_mu) + exp(log_mu)^2 * (olre.sd^2 + sd.id^2))) |>
  dplyr::select(target, model, log_mu, var) |>
  inner_join(obs) |>

  summarise(.by = model,
            c = cor(log(var), log(var.obs)))


  ggplot(aes(log(var.obs), log(var))) + geom_point() +
  facet_wrap(~ model) +

  geom_abline(slope = 1, intercept = 1)


  print()






  bind_rows(
    summaries$m1_results$evaluations |>
      mutate(model = "m1"),
    summaries$m2_results$evaluations |>
      mutate(model = "m2"),

    summaries$m4_results$evaluations |>
      mutate(model = "m4"),
    summaries$m5_results$evaluations |>
      mutate(model = "m5") ) |>

    inner_join(pois_sdid) |>

    mutate(var = if_else(is.na(olre.sd),
                         exp(log_mu) + exp(log_mu)^2 / exp(dispersion),
                         exp(log_mu) + exp(log_mu)^2 * (olre.sd^2 + sd.id^2))) |>

    mutate(dispersion = if_else(model %in% c("m1", "m2"),
                                dispersion,
                                olre.sd^2 + sd.id^2) ) |>

    dplyr::select(target, model, dispersion) |>

    pivot_wider(names_from = model, values_from = dispersion) |>


 #   summarise(c = cor(-m2, log(m5), use = "complete.obs"))

    ggplot(aes(-m2, log(m5))) + geom_point()

    dplyr::select()






estimates |>
  filter(model %in% c("m5", "m4")) |>
  filter(group == "seq_sample_id") |>
  dplyr::select(model, target, estimate) |>
  inner_join(obs) |>
  mutate(polre = mean + estimate^2) |>

  ggplot(aes(log(mean), log(polre))) + geom_point() +
  facet_wrap(~model)


  summarise(.by = model,
            c = cor(log(var), log(polre)))

  ggplot(aes(log(var), log(polre))) + geom_point() +
  facet_wrap(~model)


bind_rows(
  summaries$m1_results$evaluations |>
    mutate(model = "m1"),
  summaries$m2_results$evaluations |>
    mutate(model = "m2")) |>
  dplyr::select(model, target, dispersion) |>
  inner_join(obs) |>
  mutate(nb_var = mean + mean^2 / exp(dispersion)) |>

  ggplot(aes(log(mean), log(nb_var))) + geom_point() +
  facet_wrap(~model)


  summarise(.by = model,
            c = cor(log(var), log(nb_var)))


  ggplot(aes(log(var), log(nb_var))) + geom_point() +
  facet_wrap(~model)
