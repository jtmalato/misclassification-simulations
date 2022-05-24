#!/usr/bin/env Rscript
# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-18 10:54:09
# Title: Simulations on Steiner et al. 2020 & Cliff et al. (2019)
# Description: Using functions (simulation-functions.R) to simulate
#   1. impact of misdiagnosis in patients with infection triggered onset
# from publicly available GWAS data (Steiner et al. 2020); and
#   2. impact of misdiagnosis and misclassification in severely affected ME/CFS
# patients, considering serological data from Cliff et al. (2019).

# To run in servers -------------------------------------------------------
path <- "/tmp/RtmpvOnQRT/downloaded_packages"
dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
# Libraries ---------------------------------------------------------------
.libPaths(path)
install.packages(c("data.table", "parallel", "doParallel", "here"))
# install.packages("sendmailR", repos="http://olafmersmann.github.io/drat")
library(data.table, lib.loc = path)
library(here, lib.loc = path)
library(parallel, lib.loc = path)
library(doParallel, lib.loc = path)
registerDoParallel(20)

# Load data ---------------------------------------------------------------
steiner <- fread(here("steiner2020.csv"))
cliff <- fread(here("cliff2019.csv"))

# Load functions ----------------------------------------------------------
source(here("simulation-functions.R"))


# Run simulations ---------------------------------------------------------


# |- steiner2020 ----------------------------------------------------------

# steiner[, unique(snp)]
# steiner[cohort == "Control", n]
# steiner[cohort == "CFS", n]
# steiner[cohort == "Control", af]
# steiner[cohort == "CFS_w_ito", or]

simulation_steiner <-
  data.table(
    sim = 10000,
    n_control = steiner[cohort == "Control", allele_total],
    n_cfs = steiner[cohort == "CFS_w_ito", allele_total],
    p0 = steiner[cohort == "Control", af],
    or_t = steiner[cohort == "CFS_w_ito", or],
    se = 1,
    sp = 1
  )
simulation_steiner <- simulation_steiner[rep(seq_len(nrow(simulation_steiner)), each = length(seq(0, 1, 0.01))), ]
simulation_steiner[, gamma := rep(seq(0, 1, 0.01), length(unique(steiner$snp)))]


# |-- 1. Simulations general ----------------------------------------------

sim_steiner <-
  foreach(i = seq_len(nrow(simulation_steiner)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_steiner[i, serology_simulations(sim = sim,
                                               n_control = n_control,
                                               n_cfs = n_cfs,
                                               p0 = p0,
                                               or_t = or_t,
                                               se = se,
                                               sp = sp,
                                               gamma = gamma,
                                               test = "chisq.test")]
  }
sim_steiner_dt <- as.data.table(sim_steiner)
sim_steiner_dt[, snp := rep(unique(steiner$snp), each = length(seq(0, 1, 0.01)))]

fwrite(simulation_steiner, here("data", paste0(paste(Sys.Date(), "simulation-structure-steiner2020", sep = "_"), ".csv")))
fwrite(sim_steiner_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-steiner2020", sep = "_"), ".csv")))


# |-- 2. Simulations with fisher.test (experiments) -----------------------

sim_steiner_2 <-
  foreach(i = seq_len(nrow(simulation_steiner)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_steiner[i, serology_simulations_both_tests(sim = sim,
                                                          n_control = n_control,
                                                          n_cfs = n_cfs,
                                                          p0 = p0,
                                                          or_t = or_t,
                                                          se = se,
                                                          sp = sp,
                                                          gamma = gamma)]
  }

sim_steiner_2_dt <- as.data.table(sim_steiner_2)
sim_steiner_2_dt[, snp := rep(steiner[, as.character(unique(snp))], each = length(seq(0, 1, 0.01)))]

fwrite(sim_steiner_2_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-steiner2020-2tests", sep = "_"), ".csv")))


# |- cliff2019 ------------------------------------------------------------

# cliff[, unique(virus)]
# cliff[cohort == "Control", n]
# cliff[cohort == "CFS", n]
# cliff[cohort == "Control", pos]
# cliff[cohort == "CFS", pos]
# cliff[cohort == "CFSsa", or]

simulation_cliff <-
  data.table(
    sim = 10000,
    n_control = cliff[cohort == "Control", n],
    n_cfs = cliff[cohort == "CFS", n],
    p0 = cliff[cohort == "Control", pos],
    or_t = cliff[cohort == "CFSsa", or],
    se = 0.975,
    sp = 0.975
  )
simulation_cliff <- simulation_cliff[rep(seq_len(nrow(simulation_cliff)), each = length(seq(0, 1, 0.01))), ]
simulation_cliff[, gamma := rep(seq(0, 1, 0.01), length(unique(cliff$virus)))]


# |-- 1. Simulations general ----------------------------------------------

sim_cliff <-
  foreach(i = seq_len(nrow(simulation_cliff)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_cliff[i, serology_simulations(sim = sim,
                                             n_control = n_control,
                                             n_cfs = n_cfs,
                                             p0 = p0,
                                             or_t = or_t,
                                             se = se,
                                             sp = sp,
                                             gamma = gamma,
                                             test = "chisq.test")]
  }

sim_cliff_dt <- as.data.table(sim_cliff)
sim_cliff_dt[, virus := rep(unique(cliff$virus), each = length(seq(0, 1, 0.01)))]

fwrite(simulation_cliff, here("data", paste0(paste(Sys.Date(), "simulation-structure-cliff2019", sep = "_"), ".csv")))
fwrite(sim_cliff_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-cliff2019", sep = "_"), ".csv")))


# |-- 2. Simulations with fisher.test (experiments) -----------------------

sim_cliff_2 <-
  foreach(i = seq_len(nrow(simulation_cliff)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_cliff[i, serology_simulations_both_tests(sim = sim,
                                                        n_control = n_control,
                                                        n_cfs = n_cfs,
                                                        p0 = p0,
                                                        or_t = or_t,
                                                        se = se,
                                                        sp = sp,
                                                        gamma = gamma)]
  }

sim_cliff_2_dt <- as.data.table(sim_cliff_2)
sim_cliff_2_dt[, snp := rep(cliff[, as.character(unique(virus))], each = length(seq(0, 1, 0.01)))]

fwrite(sim_cliff_2_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-cliff2019-2tests", sep = "_"), ".csv")))

# end
