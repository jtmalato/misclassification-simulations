# Joao Malato
# 24/11/2021
# Script where simulations are produced


# Libraries ---------------------------------------------------------------
library(data.table)
library(here)


# Load functions ----------------------------------------------------------

source(here("scripts/simulation-functions.R"))


# Run simulations ---------------------------------------------------------

library(parallel)
# detectCores()
library(doParallel)
registerDoParallel(3)


# Candidate Gene ----------------------------------------------------------

simulation_gene <- CJ(
  sim    = 10000,
  n_both = c(100, 250, 500, 1000, 2500, 5000),
  p0     = c(0.05, 0.1, 0.25, 0.5),
  or_t   = c(1.25, 1.5, 2, 5, 10),
  se     = 1,
  sp     = 1,
  gamma  = seq(0, 1, 0.01)
)

sim_gene <-
  foreach(i = seq_len(nrow(simulation_gene)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_gene[i, serology_simulations(sim = sim,
                                            n_control = n_both,
                                            n_cfs = n_both,
                                            p0 = p0,
                                            or_t = or_t,
                                            se = se,
                                            sp = sp,
                                            gamma = gamma)]
  }
sim_gene_dt <- as.data.table(sim_gene)

# fwrite(simulation_gene, here("data", paste0(paste(Sys.Date(), "simulation-structure-candidate-gene", sep = "_"), ".csv")))
# fwrite(sim_gene_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-candidate-gene", sep = "_"), ".csv")))


# Serology ----------------------------------------------------------------

simulation_sero <- CJ(
  sim    = 10000,
  n_both = c(100, 250, 500, 1000, 2500, 5000),
  p0     = 0.25,
  or_t   = 1.5,
  se     = c(0.8, 0.9, 0.925, 0.975, 1.0),
  sp     = c(0.8, 0.9, 0.925, 0.975, 1.0),
  gamma  = seq(0, 1, 0.01)
)

sim_sero <-
  foreach(i = seq_len(nrow(simulation_sero)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_sero[i, serology_simulations(sim = sim,
                                                n_control = n_both,
                                                n_cfs = n_both,
                                                p0 = p0,
                                                or_t = or_t,
                                                se = se,
                                                sp = sp,
                                                gamma = gamma)]
  }
sim_sero_dt <- as.data.table(sim_sero)

# fwrite(simulation_sero, here("data", paste0(paste(Sys.Date(), "simulation-structure-serology", sep = "_"), ".csv")))
# fwrite(sim_sero_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-serology", sep = "_"), ".csv")))

# end
