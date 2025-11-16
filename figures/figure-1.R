
library(tidyverse)
library(ggtext)
library(seqwrap)
library(cowplot)


source("figures/figure-opts.R")


# Read the data from the preliminary model
m1_results <- readRDS("data/m1_results.RDS")
m1_sum <- seqwrap_summarise(m1_results, verbose = FALSE)

# Fit a trend to the dispersion data from m1
# save the data in a convenient format.
dispersion_dat <- m1_sum$evaluations


# Fit a loess model, this model can be used in the
# simulation function to represent the mean-dispersion
# relationship
trend_model_observed <- loess(dispersion ~ log_mu,
                              data = dispersion_dat,
                              span = 0.7,
                              weights = 1/(dispersion.se^2))
# Low variability scenario
trend_model_observed_low <- loess(dispersion ~ log_mu,
                              data = dispersion_dat,
                              span = 0.7)



# Plot the dispersion fit, including examples of sampling
# weights (SE of dispersion estimates)

# Select a subset for displaying se
set.seed(66)
dispersion_subset <- dispersion_dat |>
  slice_sample(n = 10)

# Dispersion model data
disp_pred_data <- data.frame(log_mu = seq(from = 0.5, to = 10, by = 0.1),
                  pred = predict(trend_model_observed,
                                 newdata = data.frame(
                                   log_mu = seq(from = 0.5,
                                                to = 10,
                                                by = 0.1))),
                  sd = trend_model_observed$s,
                  sd.low = trend_model_observed_low$s)



p1 <- dispersion_dat |>

  ggplot(aes(log_mu, dispersion)) +

  geom_point(alpha = 0.2,
             color = colors[5]) +


  geom_ribbon(data = disp_pred_data,
              aes(x = log_mu,
                  y = pred,
                  ymin = pred - sd,
                  ymax = pred + sd),
              alpha = 0.1) +

  geom_ribbon(data = disp_pred_data,
              aes(x = log_mu,
                  y = pred,
                  ymin = pred - sd.low,
                  ymax = pred + sd.low),
              alpha = 0.3) +



  geom_line(data = disp_pred_data,
            aes(log_mu, pred),
            color = colors[2],

            linewidth = 1.2) +
  # Add subset with se
 # geom_errorbar(data = dispersion_subset,
 #               aes(log_mu, dispersion,
 #                   ymin = dispersion - dispersion.se,
 #                   ymax = dispersion + dispersion.se),
 #               width = 0.2) +
#
  # Removing example data from weights calculation
 #  geom_point(data = dispersion_subset,
 #            aes(log_mu, dispersion, size = 1/dispersion.se^2),
 #            color = "black",
 #            fill = colors[5],
 #            shape = 21) +

  theme_classic() +

  labs(x = "log(&mu;)",
       y = expression(paste("log(", theta,")"))) +

  annotate("text",
           x = c(10.1, 10.1),
           y = c(3, 6.5),
           hjust = 0,
           lineheight = 0.75,
           label = c("Low variability\nscenario",
                     "High variability\nscenario"),

           size = 2.5) +

  theme(axis.title.x = element_markdown(),
        axis.title.y = element_text(),

        legend.position = "none")

p1
# Add simulated parameters

# Set the number of truly different genes
nullgenes <- 7500
condB_true <- 1250
condB_time2_true <- 1250

ngenes <- nullgenes + condB_true + condB_time2_true

  ## Set fixed (population level) effects

sim_param <- data.frame(
  beta0 = runif(nullgenes + condB_true + condB_time2_true,
                 min = 1.5,
                 max = 7),
  conditionB =  c(rep(0, nullgenes),
                   runif(
                     condB_true,
                     min = 0.2,
                     max = 1) * sample(c(-1, 1),
                                       condB_true,
                                       prob = c(0.5, 0.5),
                                       replace = TRUE),
                   rep(0, condB_time2_true)),
  timet2 = rnorm(nullgenes + condB_true + condB_time2_true,
                   0,
                   0.1),
  timet3 = rnorm(nullgenes + condB_true + condB_time2_true,
                   0,
                   0.2),
  conditionB_timet2 = rep(0, nullgenes +
                             condB_true +
                             condB_time2_true),

  conditionB_timet3 = c(rep(0, nullgenes),
                          rep(0, condB_true),
                          runif(
                            condB_true,
                            min = 0.2,
                            max = 1) * sample(c(-1, 1),
                                              condB_true,
                                              prob = c(0.25, 0.75),
                                              replace = TRUE)),
  b0_values = rlnorm(nullgenes + condB_true + condB_time2_true,
                      meanlog = -2.07,
                      sdlog = 1)) |>
  pivot_longer(names_to = "term",
               values_to = "estimate",
               cols = everything()) |>

  mutate(Term = case_when(
    term == "beta0" ~ "&beta;<sub>0</sub> (Intercept)",
    term == "conditionB" ~ "&beta;<sub>1</sub> (Group<sub>diabetes</sub>)",
    term == "timet2" ~ "&beta;<sub>2</sub> (Time<sub>2</sub>)",
    term == "timet3" ~ "&beta;<sub>3</sub> (Time<sub>3</sub>)",
    term == "conditionB_timet2" ~ "&beta;<sub>4</sub> (Time<sub>2</sub>&times;Group)",
    term == "conditionB_timet3" ~ "&beta;<sub>5</sub> (Time<sub>3</sub>&times;Group)",
    term == "b0_values" ~ "b<sub>0</sub> (SD)"
  )) |>
  mutate(type = "sim") |>
  select(Term,  estimate, type)




# Overall distributions of observed effects
p2 <- m1_sum$summaries |>
  mutate(type = "obs") |>
  select(term, estimate, type) |>
#  summarise(.by = term,
 #           s = sd(estimate)) |>
  # print()

  mutate(Term = case_when(
    term == "(Intercept)" ~ "&beta;<sub>0</sub> (Intercept)",
    term == "groupT2D" ~ "&beta;<sub>1</sub> (Group<sub>diabetes</sub>)",
    term == "timepost" ~ "&beta;<sub>2</sub> (Time<sub>2</sub>)",
    term == "timerec" ~ "&beta;<sub>3</sub> (Time<sub>3</sub>)",
    term == "timepost:groupT2D" ~ "&beta;<sub>4</sub> (Time<sub>2</sub>&times;Group)",
    term == "timerec:groupT2D" ~ "&beta;<sub>5</sub> (Time<sub>3</sub>&times;Group)",
    term == "sd__(Intercept)" ~ "b<sub>0</sub> (SD)"
  )) |>



  ggplot(aes(estimate)) +

  facet_wrap(~ Term, scales = "free", ncol = 2) +

  geom_histogram(data = sim_param,
                 aes(estimate, y = after_stat(density)),
                 alpha = 0.8,
                 boundary = 0,
                 fill = colors[2],
                 closed = "left") +

  geom_density( color = colors[5]) +

  theme_classic() +
  labs(x = "Parameter estimate") +

  theme(strip.text = element_markdown(hjust = 0),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank())


p2

# Save figure components ##################################################

figure1 <- plot_grid(plot_grid(NULL, p2,NULL,
                    ncol = 3, rel_widths = c(0.1, 1, 0)),
          plot_grid(NULL, p1, NULL,
                    rel_heights = c(0.2, 1, 0.2),
                    ncol = 1),
          ncol = 2,
          rel_widths = c(1.5,1)) +
  annotate("text", x = c(0.02,0.6),
           y = c(0.98, 0.89), label = c("A","B"))



saveRDS(figure1, "figures/figure-1.RDS")


ggsave("figures/figure-1.pdf",
       figure1,
       device = cairo_pdf,
       height = 120,
       width = 170,
       units = "mm")

