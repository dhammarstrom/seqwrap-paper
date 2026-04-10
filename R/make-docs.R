##############################################################################
#
# Make documentation
#
# This file runs files to collect data for documentation (figures, manuscript
# and supplement files). To reproduce the results presented in the manuscript
# models needs to be re-fitted. Running the full model script will
# take > 2 hours.
#
# Simulation results are re-run if make_sim is set to TRUE. Running simulation
# res takes > 24 hours (dependent on number of cores). We recommend downloading
# simulation results from Dataverse (https://doi.org/10.18710/I7U71O), see the
# Readme for instructions.
#
##############################################################################

# Re-run simulations?
make_sim <- FALSE

# Check if simulation results are present
if (!dir.exists("data_sim")) {
  warning(
    "Simulation results are not present in this repository. Results can be downloaded from Dataverse (https://doi.org/10.18710/I7U71O). The simulation results are needed to reproduce results in the manuscript."
  )
}


# Check if output folders exists
if (!dir.exists("data/")) dir.create("data/")
if (!dir.exists("data-out/")) dir.create("data-out/")

# Re-run simulations
if (make_sim) source("R/make-sims.R")

# Prepare data for models on real-world data
source("R/data-prep.R")

# Run models on real-world data
source("R/m1-m5-pillon-data.R")

# Source figure files (these are needed for the manuscript and supplement)
source("figures/figure-2.R")
source("figures/figure-3.R")
source("figures/figure-4.R")


# Render documentation
quarto::quarto_render("manuscript-docx-v3.qmd")
quarto::quarto_render("supplement.qmd")
