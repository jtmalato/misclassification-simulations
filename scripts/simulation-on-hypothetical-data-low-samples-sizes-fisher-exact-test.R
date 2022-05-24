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
  sim    = 5000,
  n_both = c(100, 250, 500),
  p0     = seq(0.05, 0.5, 0.025),
  or_t   = c(5, 10),
  se     = 1,
  sp     = 1,
  gamma  = seq(0, 1, 0.01)
)

sim_gene <-
  foreach(i = seq_len(nrow(simulation_gene)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_gene[i, serology_simulations_0(sim = sim,
                                              n_control = n_both,
                                              n_cfs = n_both,
                                              p0 = p0,
                                              or_t = or_t,
                                              se = se,
                                              sp = sp,
                                              gamma = gamma)]
  }
sim_gene_dt <- as.data.table(sim_gene)

fwrite(simulation_gene, here("data", paste0(paste(Sys.Date(), "simulation-structure-candidate-gene-chisq-fisher", sep = "_"), ".csv")))
fwrite(sim_gene_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-candidate-gene-chisq-fisher", sep = "_"), ".csv")))
