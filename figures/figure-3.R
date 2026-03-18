## Figure 3 -- Simulation results ##

library(tidyverse)
library(cowplot)
library(ggtext)
source("figures/figure-opts.R")

## Load data from simulations
source("R/simulation-functions.R")


# Extract simulations sorts out all simulations files and combines them.
# This is simulations with dispersion scenario 1.
sim_results <- extract_simulations(
  evaluations_path = "data_sim/evaluations",
  estimates_path = "data_sim/estimates",
  populationeffects_path = "data_sim/simdata/popeffect",
  disp_scenario = "s1"
)


sim_results2 <- extract_simulations(
  evaluations_path = "data_sim/evaluations2",
  estimates_path = "data_sim/estimates2",
  populationeffects_path = "data_sim/simdata2/popeffect",
  disp_scenario = "s2"
)


# Number of genes per data set with true population effect
true_effects <- bind_rows(
  sim_results$populationeffects,
  sim_results2$populationeffects
) |>
  mutate(
    effect = if_else(population_effect == 0, "true.negative", "true.positive")
  ) |>
  summarise(.by = c(term, dataset, size, disp_scenario, effect), n = n())

# Savining for manuscript
saveRDS(true_effects, "data-out/true_effects.RDS")


# Filtering for singular fits in the lmer models and pdHess/convergence in
# the other models

keep_targets <- bind_rows(sim_results$evaluations, sim_results2$evaluations) |>
  # These are not filtering away anything
  filter(!(model %in% c("m1", "m2", "m4", "m5") & pdHess == FALSE)) |>
  filter(!(model %in% c("m1", "m2", "m4", "m5") & convergence != 0)) |>

  # This has large effect on number of fitted models
  filter(!(model == "m3" & isSingular == TRUE)) |>
  dplyr::select(target, model, dataset = datasets, size, disp_scenario)


# Combine all true/false effects
est_temp <- bind_rows(sim_results$estimates, sim_results2$estimates) |>
  # Retain only terms of interest.
  filter(term %in% c("conditionB", "timet3:conditionB")) |>
  # dataset was mis-named in the simulations.
  dplyr::rename(dataset = datasets) |>
  dplyr::select(-file)


est <- est_temp |>
  inner_join(
    bind_rows(sim_results$populationeffects, sim_results2$populationeffects) |>
      mutate(target = as.character(target)) |>
      dplyr::select(-file)
  ) |>
  inner_join(
    keep_targets |>
      mutate(target = as.character(target))
  ) |>

  dplyr::select(
    target,
    term,
    estimate:p.value,
    population_effect,
    dataset,
    model,
    size,
    disp_scenario
  ) |>

  mutate(
    .by = c(model, term, size, dataset, disp_scenario),
    fdr = p.adjust(p.value, method = "fdr")
  ) |>
  mutate(
    true_effect = if_else(population_effect == 0, "neg", "pos"),
    identified_effect = if_else(fdr > 0.05, "neg", "pos"),
    true_positive = if_else(
      true_effect == "pos" &
        identified_effect == "pos",
      TRUE,
      FALSE
    ),
    false_positive = if_else(
      true_effect == "neg" &
        identified_effect == "pos",
      TRUE,
      FALSE
    )
  )


# Save est for comparison in text ###########################################

saveRDS(est, "data-out/simulation-estimates.RDS")


# Sensitivity ###############################################################

avg_sensitivity <- est |>
  filter(true_effect == "pos") |>

  summarise(
    .by = c(term, model, dataset, size, disp_scenario),
    true_positive = sum(true_positive, na.rm = TRUE)
  ) |>
  inner_join(
    true_effects |>
      filter(effect == "true.positive")
  ) |>

  mutate(sensitivity = 100 * (true_positive / n)) |>

  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),
    model_type = if_else(
      model %in% c("m1", "m2"),
      "nb",
      if_else(model %in% c("m4", "m5"), "POLRE", "normal")
    ),
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

  summarise(.by = c(model, size, term, disp_scenario), m = mean(sensitivity))


p1 <- est |>
  filter(true_effect == "pos") |>

  summarise(
    .by = c(term, model, dataset, size, disp_scenario),
    true_positive = sum(true_positive, na.rm = TRUE)
  ) |>
  inner_join(
    true_effects |>
      filter(effect == "true.positive")
  ) |>

  mutate(sensitivity = 100 * (true_positive / n)) |>

  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),
    model_type = if_else(
      model %in% c("m1", "m2"),
      "nb",
      if_else(model %in% c("m4", "m5"), "POLRE", "normal")
    ),
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

  ggplot(aes(model, sensitivity, group = paste(dataset, size))) +

  geom_line(aes(color = size), alpha = 0.7) +

  geom_point(
    data = avg_sensitivity,
    aes(model, m, group = NULL, fill = model, shape = model, alpha = size),

    color = "black",
    size = 3
  ) +

  facet_grid(term ~ disp_scenario) +

  labs(
    y = "Sensitivity (% of true effects detected)",
    fill = "Model",
    shape = "Model",
    color = "Sample size"
  ) +
  scale_fill_manual(
    values = c(colors[1], colors[2], colors[4], colors[1], colors[2])
  ) +
  scale_color_manual(values = c("gray40", "gray60", "gray80")) +
  scale_shape_manual(values = c(22, 22, 23, 25, 25)) +
  scale_alpha_manual(values = c(1, 0.75, 0.5), guide = "none") +

  theme_classic() +
  theme(
    strip.text.y = element_markdown(),
    strip.text.x = element_text(hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "white"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "left"
  )


legend <- get_plot_component(p1, "guide-box", return_all = TRUE)

p1 <- p1 + theme(legend.position = "none")

# False positives ######################################

avg_fdr <- est |>
  filter(true_effect == "neg") |>
  summarise(
    .by = c(term, model, dataset, size, disp_scenario),
    false_positive = sum(false_positive, na.rm = TRUE)
  ) |>

  inner_join(
    true_effects |>
      filter(effect == "true.negative")
  ) |>

  mutate(error = 100 * (false_positive / n)) |>

  summarise(.by = c(model, size, term, disp_scenario), m = mean(error)) |>

  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),

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


p2 <- est |>
  filter(true_effect == "neg") |>
  summarise(
    .by = c(term, model, dataset, size, disp_scenario),
    false_positive = sum(false_positive, na.rm = TRUE)
  ) |>

  inner_join(
    true_effects |>
      filter(effect == "true.negative")
  ) |>

  mutate(error = 100 * (false_positive / n)) |>

  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),

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
  ggplot(aes(model, error, group = paste(dataset, size), color = size)) +

  geom_line(aes(color = size), alpha = 0.7) +

  geom_point(
    data = avg_fdr,
    aes(model, m, group = NULL, fill = model, shape = model, alpha = size),

    color = "black",
    size = 3
  ) +

  scale_alpha_manual(values = c(1, 0.75, 0.5), guide = "none") +

  facet_grid(term ~ disp_scenario) +

  labs(
    y = "Error rate (% of true negative effects)",
    fill = "Model",
    shape = "Model",
    color = "Sample size"
  ) +
  scale_fill_manual(
    values = c(colors[1], colors[2], colors[4], colors[1], colors[2])
  ) +
  scale_color_manual(values = c("gray40", "gray60", "gray80")) +
  scale_shape_manual(values = c(22, 22, 23, 25, 25)) +

  theme_classic() +
  theme(
    strip.text.y = element_markdown(),
    strip.text.x = element_text(hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "white"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    legend.direction = "vertical"
  )


# Rank order correlation of true vs. estimate effects #####################

avg_cor <- est |>

  filter(population_effect != 0) |>

  mutate(
    .by = c(model, size, term, dataset, disp_scenario),
    rank_true = rank(population_effect),
    rank_obs = rank(estimate)
  ) |>

  summarise(
    .by = c(model, size, term, disp_scenario),
    m = cor(rank_true, rank_obs)
  ) |>
  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),

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


p3 <- est |>

  filter(population_effect != 0) |>

  mutate(
    .by = c(model, size, term, dataset, disp_scenario),
    rank_true = rank(population_effect),
    rank_obs = rank(estimate)
  ) |>

  summarise(
    .by = c(model, size, term, dataset, disp_scenario),
    cor = cor(rank_true, rank_obs)
  ) |>

  mutate(
    size = factor(
      size,
      levels = c("large", "medium", "small"),
      labels = c(
        "Large (m = 216, n = 72)",
        "Medium (m = 108, n = 36)",

        "Small (m = 54, n = 18)"
      )
    ),
    term = factor(
      term,
      levels = c("conditionB", "timet3:conditionB"),
      labels = c(
        "&beta;<sub>1</sub> group differences",
        "&beta;<sub>5</sub> group&times;time"
      )
    ),
    disp_scenario = factor(
      disp_scenario,
      levels = c("s2", "s1"),
      labels = c("Low variability", "High variability")
    ),

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

  ggplot(aes(model, cor, group = paste(dataset, size), color = size)) +

  geom_line(aes(color = size), alpha = 0.7) +
  facet_grid(term ~ disp_scenario) +

  geom_point(
    data = avg_cor,
    aes(model, m, group = NULL, fill = model, shape = model, alpha = size),

    color = "black",
    size = 3
  ) +

  scale_alpha_manual(values = c(1, 0.75, 0.5), guide = "none") +

  labs(
    y = "Correlation coefficient",
    fill = "Model",
    shape = "Model",
    color = "Sample size"
  ) +
  scale_fill_manual(
    values = c(colors[1], colors[2], colors[4], colors[1], colors[2])
  ) +
  scale_color_manual(values = c("gray40", "gray60", "gray80")) +
  scale_shape_manual(values = c(22, 22, 23, 25, 25)) +

  theme_classic() +
  theme(
    strip.text.y = element_markdown(),
    strip.text.x = element_text(hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "white"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    legend.direction = "vertical"
  )


# Combine plot ##############################################################

figure3 <- plot_grid(
  plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1)),
  plot_grid(p3, ggdraw(legend), rel_heights = c(1, 1), ncol = 1),
  rel_widths = c(1, 1)
) +
  annotate(
    "text",
    x = c(0.02, 0.02, 0.53),
    y = c(0.98, 0.49, 0.98),
    label = c("A", "C", "B")
  )


## Save RDS ##
saveRDS(figure3, "figures/figure-3.RDS")


## Save PDF ##
ggsave(
  "figures/figure-3.pdf",
  figure3,
  device = cairo_pdf,
  height = 170,
  width = 170,
  units = "mm"
)
