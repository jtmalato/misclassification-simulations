#!/usr/bin/env Rscript
# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-16 10:53:26
# Title: Functions for simulations
# Description: Script with functions used on simulations files


# Paper equations ---------------------------------------------------------

# True cases of theta1* through input of parameters theta0 OR, from equation (3)
# This parameter is dependent on 2 variables:
#   1. p0: healthy controls exposure probability
#   2. or: true odds ratio
p1_t <- function(p0, or) {
  out <- (p0 * or) / (1 - p0 + (p0 * or))
  return(out)
}
# p1_t(0.25, 1.5)
# Alternative way of coding it:
# p1_t <- function(p0, or) {
#   out <- (p0 * or) / (1 + (p0 * (or - 1)))
#   return(out)
# }

# Estimate suspected cases theta1, from equation (2)
# This parameter is dependent of 5 variables:
#   1. gamma:       misclassification rate
#   2. se:          serological test sensitivity
#   3. sp:          serological test specificity
#   4. theta0:      healthy contorls exposure probability
#   5. theta1_true: true CFS exposure probability
theta1 <- function(gamma, se, sp, theta0, theta1_true) {
  apparent_exp_pos    <- se * gamma * theta0
  apparent_nonexp_pos <- (1-sp) * gamma * (1-theta0)
  true_exp_pos        <-  se * (1-gamma) * theta1_true
  tru_nonexp_pos      <- (1-sp) * (1-gamma) * (1-theta1_true)
  return(apparent_exp_pos + apparent_nonexp_pos + true_exp_pos + tru_nonexp_pos)
}

# Order and estimate p-values form the 2-by-2 contingency table observed
# Using the chi-squared test
get_p_value_chisq <- function(data) {
  return(suppressWarnings(chisq.test(matrix(data, ncol = 2), correct = FALSE)$p.value))
}
# Order and estimate p-values form the 2-by-2 contingency table observed
# Using the Fisher's exact test
get_p_value_fisher <- function(data) {
  return(suppressWarnings(fisher.test(matrix(data, ncol = 2))$p.value))
}


# Estimate confidence interval based on vector of p-values
confidence_interval <- function(vector, interval = 0.95) {
  # standard deviation of sample
  v_sd <- sd(vector)
  # sample size
  v_n <- length(vector)
  # mean of sample
  v_mean <- mean(vector)
  # error according to t distribution
  v_error <- qt((interval + 1)/2, df = v_n - 1) * v_sd / sqrt(v_n)
  # confidence interval as a vector
  return(c(v_mean - v_error, v_mean + v_error))
}


# Simulation structure ----------------------------------------------------

# Simulations using using only chisq.test or fisher.test (returning p-value CIs)
serology_simulations <- function(
  sim       = 10000,
  n_control = 100,
  n_cfs     = 100,
  p0        = 0.5,
  or_t      = 2,
  se        = 1,
  sp        = 1,
  gamma     = 0.5,
  test      = "chisq.test"
)
{

  if(test != "chisq.test" || test != "fisher.test") {
    warning(
      paste0("Wrong test used: '", test, "'. Must choose between 'chisq.test' or 'fisher.test'."),
      "\n",
      "'chisq.test' used to perform simulations."
    )
    test <- "chisq.test"
  }

  require(data.table)
  p1_true <- p1_t(p0, or_t)

  # Controls ------------------------------------------
  # true exposed/non-exposed healthy controls
  xp_c <- data.table::data.table(xp = rbinom(sim, n_control, p0))[, nxp := n_control-xp][]
  # structure of serological test
  outcome_c <- xp_c[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_c$xp-t_pos,
                                  f_pos = xp_c$nxp-t_neg)][]
  # serological test outcome
  serotest_c <- outcome_c[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]

  # CFS -----------------------------------------------
  # apparent/false
  n_apparent <- rbinom(sim, n_cfs, gamma)
  xp_f <- data.table::data.table(xp = rbinom(sim, n_apparent, p0))[, nxp := n_apparent-xp][]
  # structure of serological test
  outcome_f <- xp_f[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_f$xp-t_pos,
                                  f_pos = xp_f$nxp-t_neg)][]
  # serological test outcome
  serotest_f <- outcome_f[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]
  # true
  n_true <- n_cfs - n_apparent
  xp_t <- data.table::data.table(xp = rbinom(sim, n_true, p1_true))[, nxp := n_true-xp][]
  # structure of serological test
  outcome_t <- xp_t[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_t$xp-t_pos,
                                  f_pos = xp_t$nxp-t_neg)][]
  # serological test outcome
  serotest_t <- outcome_t[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]

  # Observed ------------------------------------------
  # gather data created
  data_obs <- cbind(hc = serotest_c, cfs = serotest_f + serotest_t)

  if(test == "chisq.test") {
    vector_p_value <- apply(data_obs, 1, get_p_value_chisq)
    p_value <- vector_p_value < 0.05
    if(sum(is.na(p_value)) != 0) {p_value[is.na(p_value)] <- FALSE}
    ci <- confidence_interval(p_value)
  } else { # test = "fisher.test"
    vector_p_value <- apply(data_obs, 1, get_p_value_fisher)
    p_value <- vector_p_value < 0.05
    if(sum(is.na(p_value)) != 0) {p_value[is.na(p_value)] <- FALSE}
    ci <- confidence_interval(p_value)
  }

  # overall exposure probability on any suspected case
  p1 <- theta1(gamma = gamma, se = se, sp = sp, theta0 = p0, theta1_true = p1_true)
  # estimated OR of cases against hc based on risk factor exposure
  odds <- (p1 * (1-p0)) / ((1-p1) * p0)

  out <- c(
    sim         = sim,
    sample_size = n_control,
    sensitivity = se,
    specificity = sp,
    misrate     = gamma,
    p           = mean(p_value),
    p_min       = ci[1],
    p_max       = ci[2],
    overall_or  = odds,
    or_t        = or_t,
    theta0      = p0,
    theta1_true = p1_true,
    theta1      = p1
  )
  return(out)
}


# function to run simulations using both chisq.test and fisher.test (no CIs)
serology_simulations_both_tests <- function(
  sim       = 10000,
  n_control = 100,
  n_cfs     = 100,
  p0        = 0.5,
  or_t      = 2,
  se        = 1,
  sp        = 1,
  gamma     = 0.5
  )
  {
  require(data.table)
  p1_true <- p1_t(p0, or_t)

  # Controls ------------------------------------------
  # true exposed/non-exposed healthy controls
  xp_c <- data.table::data.table(xp = rbinom(sim, n_control, p0))[, nxp := n_control-xp][]
  # structure of serological test
  outcome_c <- xp_c[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_c$xp-t_pos,
                                  f_pos = xp_c$nxp-t_neg)][]
  # serological test outcome
  serotest_c <- outcome_c[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]

  # CFS -----------------------------------------------
  # apparent/false
  n_apparent <- rbinom(sim, n_cfs, gamma)
  xp_f <- data.table::data.table(xp = rbinom(sim, n_apparent, p0))[, nxp := n_apparent-xp][]
  # structure of serological test
  outcome_f <- xp_f[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_f$xp-t_pos,
                                  f_pos = xp_f$nxp-t_neg)][]
  # serological test outcome
  serotest_f <- outcome_f[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]

  # true
  n_true <- n_cfs - n_apparent
  xp_t <- data.table::data.table(xp = rbinom(sim, n_true, p1_true))[, nxp := n_true-xp][]
  # structure of serological test
  outcome_t <- xp_t[, .(t_pos = rbinom(sim, xp, se),
                        t_neg = rbinom(sim, nxp, sp))][
                          , `:=` (f_neg = xp_t$xp-t_pos,
                                  f_pos = xp_t$nxp-t_neg)][]
  # serological test outcome
  serotest_t <- outcome_t[, .(pos = t_pos + f_pos,
                              neg = t_neg + f_neg)]

  # Observed ------------------------------------------
  # gather data created
  data_obs <- cbind(hc = serotest_c, cfs = serotest_f + serotest_t)

  # test H0
  vector_p_value_chisq <- apply(data_obs, 1, get_p_value_chisq)
  vector_p_value_fisher <- apply(data_obs, 1, get_p_value_fisher)
  # power to reject H0
  p_value_chisq <- vector_p_value_chisq < 0.05
  p_value_fisher <- vector_p_value_fisher < 0.05
  if(sum(is.na(p_value_chisq)) != 0) {p_value_chisq[is.na(p_value_chisq)] <- FALSE}
  if(sum(is.na(p_value_fisher)) != 0) {p_value_fisher[is.na(p_value_fisher)] <- FALSE}
  # ci <- confidence_interval(p_value)

  # overall exposure probability on any suspected case
  p1 <- theta1(gamma=gamma, se=se, sp=sp, theta0=p0, theta1_true=p1_true)

  # estimated OR of cases agains hc based on risk factor exposure
  odds <- (p1 * (1-p0)) / ((1-p1) * p0)

  # out
  out <- c(
    sim         = sim,
    sample_size = n_control,
    sensitivity = se,
    specificity = sp,
    misrate     = gamma,
    p_chisq     = mean(p_value_chisq),
    p_fisher    = mean(p_value_fisher),
    # p_min       = ci[1],
    # p_max       = ci[2],
    overall_or  = odds,
    or_t        = or_t,
    theta0      = p0,
    theta1_true = p1_true,
    theta1      = p1
  )
  return(out)
}

# end
