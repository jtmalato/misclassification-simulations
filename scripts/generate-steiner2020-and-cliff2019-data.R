#!/usr/bin/env Rscript
# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-18 10:56:29
# Title: Collect data from Steiner et al. 2020 & Cliff et al. (2019)
# Description: Using published data from Steiner et al. 2020 and Cliff et al. (2019),
# I collected available values from Tables 3 and 1, respectively.


# Libraries ---------------------------------------------------------------
library(here)
library(data.table)


# Steiner 2020 ------------------------------------------------------------


# |- Generate data --------------------------------------------------------

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

# estimate odds ratio sandard error
for(i in seq_along(unique(steiner$cohort))) {
  for(j in unique(steiner$snp)) {
    steiner[snp == j & cohort %in% c(unique(steiner$cohort)[c(i, 4)]), or_se := sqrt(sum(1/c(allele_normal, allele_risk)))]
  }
}

# odds ratio 95% confidence interval
steiner[, `:=` (or_lower = exp(c(log(or) - qnorm(1 - 0.05/2) * or_se)),
                or_upper = exp(c(log(or) + qnorm(1 - 0.05/2) * or_se)))]

# Two-tailed test p-value based on OR
steiner[, or_p := 2 * pnorm(abs(log(or)) / or_se, lower.tail = FALSE)]
# steiner[, pnorm(abs(log(or)) / or_se, lower.tail = TRUE)]
# steiner[, pnorm(abs(log(or)) / or_se, lower.tail = FALSE)]
# steiner[, 2 * pnorm(abs(log(or)) / or_se, lower.tail = FALSE)]

# for(i in seq_along(unique(steiner$cohort))) {
#   for(j in unique(steiner$snp)) {
#     ps <- chisq.test(steiner[snp == j & cohort %in% c(unique(steiner$cohort)[c(i, 4)]), matrix(c(allele_risk, allele_normal), ncol = 2, byrow = TRUE)])$p.value
#     steiner[snp == j & cohort %in% c(unique(steiner$cohort)[c(i, 4)]), chisq_p := ps]
#   }
# }
# steiner[cohort == "Control", chisq_p := 1]


# |- Using patients with infection triggered onset ------------------------

# (CFS_w_ito)

steiner_tab_cfs <- steiner[cohort %in% c("CFS_w_ito"), .(snp, af = round(af, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3), order = seq_along(unique(steiner$snp)))]
steiner_tab_hc <- steiner[cohort %in% c("Control"), .(snp, af = round(af, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3))]

steiner_tab_cfs[order(or_p), -"order"]
steiner_tab_hc[steiner_tab_cfs[order(or_p)]$order]


# Cliff 2019 --------------------------------------------------------------


# |- Generate data --------------------------------------------------------

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
# odds/likelihood of pos
cliff[, odds := pos / neg]
# odds ratio comparing each ME/CFS category against healthy controls
cliff[, or := odds / odds[4], by = virus]

# # Alternative: use Rfast package to estimate OR, Cis and p-values
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

# for(i in seq_along(unique(cliff$cohort))) {
#   for(j in unique(cliff$virus)) {
#     ps <- chisq.test(cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]), matrix(c(exposed, nexposed), ncol = 2, byrow = TRUE)])$p.value
#     cliff[virus == j & cohort %in% c(unique(cliff$cohort)[c(i, 4)]), chisq_p := ps]
#   }
# }
# cliff[cohort == "Control", chisq_p := 1]


# |-- Using severely affected patients ------------------------------------

# (CFSsa)

cliff_tab_cfs <- cliff[cohort %in% c("CFSsa"), .(virus, pos = round(pos, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3), order = seq_along(unique(cliff$virus)))]
cliff_tab_hc <- cliff[cohort %in% c("Control"), .(virus, pos = round(pos, 2), or = round(or, 2), or_lower = round(or_lower, 2), or_upper = round(or_upper, 2), or_p = round(or_p, 3))]

cliff_tab_cfs[order(or_p), -"order"]
cliff_tab_hc[cliff_tab_cfs[order(or_p)]$order]


# Save data ---------------------------------------------------------------

# save steiner2020 table for simulations
fwrite(steiner, here("data", "steiner2020.csv"))
# save cliff2019 table for simulations
fwrite(cliff, here("data", "cliff2019.csv"))

# end
