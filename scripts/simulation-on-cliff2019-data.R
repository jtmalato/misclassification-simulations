# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-16 11:50:02
# Title: Simluations using Cliff et al. (2019) data
# Description: Script where simulations on serology study Cliff et al. (2019) are produced


# Libraries ---------------------------------------------------------------
library(here)
library(data.table)


# Load functions ----------------------------------------------------------

source(here("scripts/simulation-functions.R"))


# Real world simulation ---------------------------------------------------


# Cliff 2019 --------------------------------------------------------------

## Estimate table values --------------------------------------------------

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
# number of non-exposed individuals
cliff[, nexposed := n - exposed]
# proportions of positive and negative individuals
cliff[, pos := exposed / n][, neg := 1 - pos]
# odds (i.e., likelihood) of being positive
cliff[, odds := pos / neg]
# odds ratio comparing ME/CFS against healthy controls
cliff[, or := odds / odds[4], by = virus]

# # use Rfast package to estimate OR, Cis and p-values
# for(i in seq_along(unique(cliff$cohort))) {
#   for(j in unique(cliff$virus)) {
#     cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]),
#           `:=` (
#             or_pckg = Rfast::odds.ratio(matrix(c(exposed, nexposed), ncol = 2, byrow = FALSE))[[1]][1],
#             or_lower = Rfast::odds.ratio(matrix(c(exposed, nexposed), ncol = 2, byrow = FALSE))[[2]][1],
#             or_upper = Rfast::odds.ratio(matrix(c(exposed, nexposed), ncol = 2, byrow = FALSE))[[2]][2],
#             or_p = Rfast::odds.ratio(matrix(c(exposed, nexposed), ncol = 2, byrow = FALSE))[[1]][2])
#           ]
#   }
# }

# odds ratio standard error
for(i in seq_along(unique(cliff$cohort))) {
  for(j in unique(cliff$virus)) {
    cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]), or_se := sqrt(sum(1/c(exposed, nexposed)))]
  }
}
# odds ratio 95% confidence interval
cliff[, `:=` (or_lower = exp(log(or) - qnorm(1 - 0.05/2) * or_se),
              or_upper = exp(log(or) + qnorm(1 - 0.05/2) * or_se))]
# estimate contigency table (with healthy control) p-value
cliff[, or_p := 2 * pnorm(abs(log(or)) / or_se, lower.tail = FALSE)]

for(i in seq_along(unique(cliff$cohort))) {
  for(j in unique(cliff$virus)) {
    ps <- chisq.test(cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]), matrix(c(exposed, nexposed), ncol = 2, byrow = TRUE)])$p.value
    cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]), chisq_p := ps]
  }
}
cliff[cohort == "Control", chisq_p := 1]



cliff[cohort %in% c("CFSsa"), .(virus, pos = round(pos, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3), chisq_p = round(chisq_p, 3))][c(3, 4, 2, 1, 5, 6)]
cliff[cohort %in% c("Control"), .(virus, pos = round(pos, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3), chisq_p = round(chisq_p, 3))][c(3, 4, 2, 1, 5, 6)]


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
                                             gamma = gamma,
                                             test = "chisq.test")]
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
