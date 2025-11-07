
## Figure 2 -- Simulation results ##


library(tidyverse)


## Load data from simulations
evaluations <- readRDS("data/sim_evaluations.RDS")
effects <-  readRDS("data/sim_populationeffects.RDS")


estimate_files <- list.files("data/estimates/")
estimates <- list()
for(i in seq_along(estimate_files))  {

  estimates[[i]] <- readRDS(paste0("data/estimates/", estimate_files[i]))

}

estimates <- bind_rows(estimates)




temp <- estimates |>
  filter(term %in% c("conditionB", "timet3:conditionB")) |>
  inner_join(effects |>
               mutate(target = as.character(target))) |>

  select(target, term, estimate:p.value,population_effect, dataset, model, size) |>

  mutate(.by = c(model, term, size, dataset),
         fdr = p.adjust(p.value, method = "fdr")) |>
  mutate(true_effect = if_else(population_effect == 0, "neg", "pos"),
         identified_effect = if_else(fdr > 0.05, "neg", "pos"),
         true_positive = if_else(true_effect == "pos" &
                                   identified_effect == "pos", TRUE, FALSE),
         false_positive = if_else(true_effect == "neg" &
                                    identified_effect == "pos", TRUE, FALSE)) |>

  print()
