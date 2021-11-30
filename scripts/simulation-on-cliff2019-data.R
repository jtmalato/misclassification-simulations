# Joao Malato
# 29/11/2021
# Script where simulations on serology study Cliff et al. () are produced


# Libraries ---------------------------------------------------------------
library(data.table)
library(here)

# Load functions ----------------------------------------------------------

source(here("scripts/simulation-functions.R"))


# real world simulation ---------------------------------------------------


# Cliff 2019 --------------------------------------------------------------


# work with the significant ones
cliff <- data.table(virus = rep(c("CMV", "EBV", "HSV1", "HSV2", "VZV", "HHV6"), each = 4),
                    cohort = rep(c("CFSsa", "CFSmm", "CFS", "Control"), 6),
                    n = rep(c(54, 197, 54+197, 107), 6),
                    exposed = c(18, 57, 18+57, 40,
                                48, 171, 48+171, 99,
                                29, 82, 29+82, 45,
                                22, 76, 22+76, 36,
                                52, 190, 52+190, 104,
                                52, 177, 52+177, 102))
# dcast(cliff, cohort + n ~ virus, value.var = "exposed")
# cliff <- cliff[!cohort %in% c("CFSmm", "CFSsa")]
cliff[, pos := exposed / n][, neg := 1 - pos]
cliff[, odds := pos / neg]
cliff[, or := odds / odds[4], by = virus]


cliff[, virus := factor(virus, levels = c("CMV", "EBV", "HSV1", "HSV2", "VZV", "HHV6"))]
dcast(cliff[cohort %in% c("CFS", "Control"), .(virus, cohort, pos = round(pos, 2), n, or = round(or, 2))], virus ~ cohort, value.var = "pos")
dcast(cliff[cohort %in% c("CFS", "Control"), .(virus, cohort, pos = round(pos, 2), n, or = round(or, 2))], virus ~ cohort, value.var = "or")
dcast(cliff[cohort %in% c("CFSsa", "Control"), .(virus, cohort, pos = round(pos, 2), n, or = round(or, 2))], virus ~ cohort, value.var = "or")
# dcast(cliff[cohort %in% c("CFS", "Control"), .(snp, cohort, af = round(af, 2), n, or = round(or, 2))], snp ~ cohort, value.var = "af")
# round(cbind(t1=cliff[cohort %in% c("CFS"), af],
#       t1_estimate = p1_t(cliff[cohort %in% c("Control"), af], cliff[cohort %in% c("CFS_w_ito"), or])), 2)
# Run simulations ---------------------------------------------------------

cliff[, unique(virus)]
cliff[cohort == "Control", n]
cliff[cohort == "CFS", n]
cliff[cohort == "Control", pos]
cliff[cohort == "CFS", pos]
cliff[cohort == "CFSsa", or]

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

library(parallel)
# detectCores()
library(doParallel)
registerDoParallel(3)

sim_cliff <-
  foreach(i = seq_len(nrow(simulation_cliff)), .combine = rbind, .packages = "data.table") %dopar% {
    simulation_cliff[i, serology_simulations(sim = sim,
                                               n_control = n_control,
                                               n_cfs = n_cfs,
                                               p0 = p0,
                                               or_t = or_t,
                                               se = se,
                                               sp = sp,
                                               gamma = gamma)]
  }

sim_cliff_dt <- as.data.table(sim_cliff)
sim_cliff_dt[, virus := rep(unique(cliff$virus), each = length(seq(0, 1, 0.01)))]

# fwrite(simulation_cliff, here("data", paste0(paste(Sys.Date(), "simulation-structure-cliff2019", sep = "_"), ".csv")))
# fwrite(sim_cliff_dt, here("data", paste0(paste(Sys.Date(), "simulation-results-cliff2019", sep = "_"), ".csv")))
# end

pnorm(1, 0, 1) - pnorm(-1, 0, 1)
pnorm(2, 0, 1) - pnorm(-2, 0, 1)
pnorm(3, 0, 1) - pnorm(-3, 0, 1)
pnorm(4, 0, 1) - pnorm(-4, 0, 1)
