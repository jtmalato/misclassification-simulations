# Joao Malato
# 29/11/2021
# Script where simulations on GWAS Steiner et al. 2020 are produced


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


# Steiner 2020 ------------------------------------------------------------

# steiner table results
steiner <-
  data.table(
    snp = rep(c("PTPN22", "CTLA4", "IRF5", "TNF1", "TNF2"), each = 4),
    cohort = rep(c("CFS", "CFS_w_ito", "CFS_wo_ito", "Control"), 5),
    homoz_risk = c(2, 1, 1, 2, 114, 97, 17, 61, 74, 63, 11, 50, 6, 5, 1, 4, 5, 5, 0, 2),
    heteroz = c(68, 57, 11, 29, 155, 112, 43, 102, 140, 102, 38, 104, 77, 56, 21, 55, 52, 41, 11, 47),
    homoz_normal = c(235, 174, 61, 170, 36, 23, 13, 38, 91, 67, 24, 47, 222, 171, 51, 142, 248, 186, 62, 150)
  )

# add variables
steiner[, `:=` (
  allele_risk = (homoz_risk * 2) + heteroz,
  allele_normal = (homoz_normal * 2) + heteroz)
  ][, allele_total := allele_risk + allele_normal
    ][, `:=` (af = allele_risk / allele_total, n = allele_total / 2)
      ][, odds := af / (1 - af)
        ][, or := odds / odds[4], by = snp]
# fwrite(steiner, here("data", "candidate-gene-steiner2020.csv"))



steiner[, snp := factor(snp, levels = c("PTPN22", "CTLA4", "TNF1", "TNF2", "IRF5"))]
# dcast(steiner[cohort %in% c("CFS", "Control"), .(snp, cohort, af = round(af, 2), n, or = round(or, 2))], snp ~ cohort, value.var = "af")
# dcast(steiner[cohort %in% c("CFS", "Control"), .(snp, cohort, af = round(af, 2), n, or = round(or, 2))], snp ~ cohort, value.var = "or")
dcast(steiner[cohort %in% c("CFS_w_ito", "Control"), .(snp, cohort, af = round(af, 2), n, or = round(or, 2))], snp ~ cohort, value.var = "or")
# dcast(steiner[cohort %in% c("CFS", "Control"), .(snp, cohort, af = round(af, 2), n, or = round(or, 2))], snp ~ cohort, value.var = "af")
# round(cbind(t1=steiner[cohort %in% c("CFS"), af],
#       t1_estimate = p1_t(steiner[cohort %in% c("Control"), af], steiner[cohort %in% c("CFS_w_ito"), or])), 2)


# Run simulations ---------------------------------------------------------

library(parallel)
# detectCores()
library(doParallel)
registerDoParallel(3)


# 1) ----------------------------------------------------------------------

steiner[, unique(snp)]
steiner[cohort == "Control", n]
steiner[cohort == "CFS", n]
steiner[cohort == "Control", af]
steiner[cohort == "CFS_w_ito", or]

simulation_steiner <-
  data.table(
    sim = 10000,
    n_control = steiner[cohort == "Control", n],
    n_cfs = steiner[cohort == "CFS", n],
    p0 = steiner[cohort == "Control", af],
    or_t = steiner[cohort == "CFS_w_ito", or],
    se = 1,
    sp = 1
  )

simulation_steiner <- simulation_steiner[rep(seq_len(nrow(simulation_steiner)), each = length(seq(0, 1, 0.01))), ]
simulation_steiner[, gamma := rep(seq(0, 1, 0.01), length(unique(steiner$snp)))]

sim_steiner <-
  foreach(i = seq_len(nrow(simulation_steiner)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_steiner[i, serology_simulations(sim = sim,
                                               n_control = n_control,
                                               n_cfs = n_cfs,
                                               p0 = p0,
                                               or_t = or_t,
                                               se = se,
                                               sp = sp,
                                               gamma = gamma)]
  }
sim_steiner_dt <- as.data.table(sim_steiner)
sim_steiner_dt[, snp := rep(unique(steiner$snp), each = length(seq(0, 1, 0.01)))]

# fwrite(simulation_steiner, here("data", paste0(paste(Sys.Date(), "simulation-structure-steiner2020", sep = "_"), ".csv")))
# fwrite(sim_steiner_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-steiner2020", sep = "_"), ".csv")))


# 2)  ---------------------------------------------------------------------

steiner[, levels(snp)]
steiner[cohort == "Control", n]
steiner[cohort == "CFS", n]
steiner[cohort == "Control", af]
steiner[cohort == "CFS_wo_ito", af]
steiner[cohort == "CFS_w_ito", or]


simulation_steiner <-
  data.table(
    sim = 10000,
    n_control = steiner[cohort == "Control", n],
    n_cfs = steiner[cohort == "CFS", n],
    p0 = steiner[cohort == "Control", af],
    p0_cfs = steiner[cohort == "CFS_wo_ito", af],
    or_t = steiner[cohort == "CFS_w_ito", or],
    se = 1,
    sp = 1
  )

simulation_steiner <- simulation_steiner[rep(seq_len(nrow(simulation_steiner)), each = length(seq(0, 1, 0.01))), ]
simulation_steiner[, gamma := rep(seq(0, 1, 0.01), length(levels(steiner$snp)))]

sim_steiner <-
  foreach(i = seq_len(nrow(simulation_steiner)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_steiner[i, serology_simulations_2(sim = sim,
                                                 n_control = n_control,
                                                 n_cfs = n_cfs,
                                                 p0 = p0,
                                                 p0_cfs = p0_cfs,
                                                 or_t = or_t,
                                                 se = se,
                                                 sp = sp,
                                                 gamma = gamma)]
  }
sim_steiner_dt <- as.data.table(sim_steiner)
sim_steiner_dt[, snp := rep(steiner[, as.character(unique(snp))], each = length(seq(0, 1, 0.01)))]

# fwrite(simulation_steiner, here("data", paste0(paste(Sys.Date(), "simulation-structure-steiner2020-2", sep = "_"), ".csv")))
# fwrite(sim_steiner_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-steiner2020-2", sep = "_"), ".csv")))
