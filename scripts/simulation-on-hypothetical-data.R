#!/usr/bin/env Rscript
# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-16 10:51:37
# Title: simulations on hypothetical data for misclassification paper.
# Description: Make use of the `simulation-functions.R` document and generate
# simulations on the two hypothetical scenarios described in the paper.

# To run in servers -------------------------------------------------------
path <- "/tmp/RtmpvOnQRT/downloaded_packages"
dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
# Libraries ---------------------------------------------------------------
.libPaths(path)
install.packages(c("data.table", "parallel", "doParallel", "here"))

# library(here)
library(data.table, lib.loc = path)
library(here, lib.loc = path)
library(parallel, lib.loc = path)
library(doParallel, lib.loc = path)
registerDoParallel(20)

# Load functions ----------------------------------------------------------
source(here("simulation-functions.R"))


# Run simulations ---------------------------------------------------------


# |- Candidate Gene -------------------------------------------------------

simulation_gene <- CJ(
  sim    = 10000,
  n_both = c(100, 250, 500, 1000, 2500, 5000),
  p0     = c(0.01, 0.05, 0.1, 0.25, 0.5), # last time we tried to use p0 = 0.01
  or_t   = c(1.25, 1.5, 2, 3, 5, 10),
  se     = 1,
  sp     = 1,
  gamma  = seq(0, 1, 0.01)
)

sim_gene <-
  foreach(i = seq_len(nrow(simulation_gene)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_gene[i, serology_simulations(sim       = sim,
                                            n_control = n_both,
                                            n_cfs     = n_both,
                                            p0        = p0,
                                            or_t      = or_t,
                                            se        = se,
                                            sp        = sp,
                                            gamma     = gamma,
                                            test      = "chisq.test")]
  }
sim_gene_dt <- as.data.table(sim_gene)

fwrite(simulation_gene, here("data", paste0(paste(Sys.Date(), "simulation-structure-candidate-gene", sep = "_"), ".csv")))
fwrite(sim_gene_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-candidate-gene", sep = "_"), ".csv")))


# |- Serology -------------------------------------------------------------

simulation_sero <- CJ(
  sim    = 10000,
  n_both = c(100, 250, 500, 1000, 2500, 5000),
  # p0     = c(0.01, 0.05, 0.1, 0.25, 0.5),
  p0     = 0.25,
  # or_t   = c(1.25, 1.5, 2, 3, 5, 10),
  or_t   = 3,
  se     = c(0.8, 0.9, 0.925, 0.975, 1.0),
  sp     = c(0.8, 0.9, 0.925, 0.975, 1.0),
  gamma  = seq(0, 1, 0.01)
)

sim_sero <-
  foreach(i = seq_len(nrow(simulation_sero)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_sero[i, serology_simulations(sim       = sim,
                                            n_control = n_both,
                                            n_cfs     = n_both,
                                            p0        = p0,
                                            or_t      = or_t,
                                            se        = se,
                                            sp        = sp,
                                            gamma     = gamma,
                                            test      = "chisq.test")]
  }
sim_sero_dt <- as.data.table(sim_sero)

fwrite(simulation_sero, here("data", paste0(paste(Sys.Date(), "simulation-structure-serology-or", sep = "_"), ".csv")))
fwrite(sim_sero_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-serology-or", sep = "_"), ".csv")))

# end
