#!/usr/bin/env Rscript
# Header ------------------------------------------------------------------
# Author: Joao Malato
# Date: 2022-05-18 12:41:15
# Title: Analysis of simulation results
# Description: Script for figures and tables


# Libraries ---------------------------------------------------------------
library(extrafont)
loadfonts(device = "win")
library(data.table)
library(here)
library(magrittr)
library(ggplot2)
library(lemon)
library(viridis)
library(MetBrewer)
library(colorspace)


# Load data ---------------------------------------------------------------

# simulation - gene scenario
dt1 <- fread(here("data-cluster/no-correction/2022-05-18_simulation-results-candidate-gene.csv"))
# dt1[, lapply(.SD, \(x) length(unique(x)))]
dt1[, ss_f := factor(sample_size)]
dt1[, ss_f := forcats::fct_rev(ss_f)]
dt1[, or_f := factor(or_t, levels = rev(as.character(unique(or_t))))]
# don't consider theta0 = 0.01
dt1 <- dt1[theta0 != 0.01]
dt1[, theta0_f := factor(theta0, levels = c(0.5, 0.25, 0.1, 0.05))]

# simulation - serological scenario
dt2 <- fread(here("data-cluster/no-correction/2022-05-18_simulation-results-serology-or.csv"))
dt2[, ss_f := factor(sample_size)]
dt2[, `:=` (se_f = factor(sensitivity, levels = rev(as.character(unique(sensitivity)))),
            sp_f = factor(specificity, levels = rev(as.character(unique(specificity)))))]


# Theme -------------------------------------------------------------------

theme_png <- function() {
  theme_classic(base_size = 18, base_family = "TT Arial") +
    theme(
      text = element_text(size = 18,  family =),
      panel.border = element_blank(),
      # panel.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted", size = 0.4),
      panel.grid.major.y = element_line(colour = c("gray", NA, rep("gray", 5)), linetype = "dotted", size = 0.4),
      axis.title.x = element_text(vjust = -2.5),
      axis.title.y = element_text(vjust = 4),
      strip.background = element_blank(),
      strip.text = element_text(size = 20),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.background = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.spacing.x = unit(3, 'pt'),
      legend.margin = margin(0, 0, -4, 0),
      panel.spacing.x = unit(1, "cm"),
      panel.spacing.y = unit(1, "cm"),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
}

theme_pdf <- function() {
  theme_classic(base_size = 18, base_family = "LM Roman 10") +
    theme(
      text = element_text(size = 18,  family =),
      panel.border = element_blank(),
      # panel.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted", size = 0.4),
      panel.grid.major.y = element_line(colour = c("gray", NA, rep("gray", 5)), linetype = "dotted", size = 0.4),
      axis.title.x = element_text(vjust = -2.5),
      axis.title.y = element_text(vjust = 4),
      strip.background = element_blank(),
      strip.text = element_text(size = 20),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.background = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.spacing.x = unit(3, 'pt'),
      legend.margin = margin(0, 0, -4, 0),
      panel.spacing.x = unit(1, "cm"),
      panel.spacing.y = unit(1, "cm"),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
}

theme_png_no_background <- function() {
  theme_classic(base_size = 18, base_family = "TT Arial") +
    theme(
      text = element_text(size = 18,  family =),
      panel.border = element_blank(),
      plot.background = element_blank(),
      # panel.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted", size = 0.4),
      panel.grid.major.y = element_line(colour = c("gray", NA, rep("gray", 5)), linetype = "dotted", size = 0.4),
      axis.title.x = element_text(vjust = -2.5),
      axis.title.y = element_text(vjust = 4),
      strip.background = element_blank(),
      strip.text = element_text(size = 20),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.background = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.spacing.x = unit(3, 'pt'),
      legend.margin = margin(0, 0, -4, 0),
      panel.spacing.x = unit(1, "cm"),
      panel.spacing.y = unit(1, "cm"),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
}


# Functions ---------------------------------------------------------------

gg_base <- function(dt, if_p = TRUE) {
  if(if_p) {
    gg_aes <- ggplot(dt, aes(misrate, p))
  } else {
    gg_aes <- ggplot(dt, aes(misrate))
  }
  gg_aes +
    # annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
    scale_y_continuous(limits = c(0,1),
                       breaks = sort(c(0.05, seq(0,1, 0.2))),
                       labels = c("", expression(alpha), "0.2","0.4","0.6",expression(1-beta), "1")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1")) +
    coord_cartesian(expand = FALSE, clip = "off") +
    labs(x = expression(paste("Misclassification rate,")~gamma),
         y = "Probability of rejecting null hypothesis",
         colour = expression(paste("Sample size")~(n[0]==n[1])))
}

prepare_ribbon_data <- function(dt,
                                sample_sizes = 1000,
                                sensitivitys = 1,
                                specificitys = 1,
                                or_ts = 5,
                                theta0s = 0.5,
                                bys = "sample_size") {
  dt_list <-
    dt[sample_size %in% sample_sizes &
         sensitivity %in% sensitivitys &
         specificity %in% specificitys &
         or_t %in% or_ts &
         theta0 %in% theta0s] %>%
    split(by = bys)

  dt_out <- merge(dt_list[[1]][, .(sample_size, misrate, p, or_t, theta0)],
                  dt_list[[2]][, .(sample_size, misrate, p, or_t, theta0)],
                  by = "misrate",
                  suffixes = paste0("_", names(dt_list)))
  return(dt_out)
}


# Figure parameters (to save) ---------------------------------------------

scale <- 1
dpi <- 320
width <- 13
height <- 16


# Simulation results ------------------------------------------------------


# |- Candidate gene -------------------------------------------------------


# |-- 1. Figures ----------------------------------------------------------

gg1 <-
  dt1 %>%
  ggplot(aes(misrate, p)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6",expression(1-beta), "1")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis",
       colour = expression(paste("Sample size")~(n[0]==n[1])~paste("  "))) +
  facet_grid(or_f ~ theta0_f,
             labeller = label_bquote(
               rows = Delta[T] == .(as.character(or_f)),
               cols = theta[0] == .(as.character(theta0_f)))
  ) +
  geom_line(aes(col = ss_f), lwd = 1.2) +
  # scale_colour_manual(values = darken(viridis(6), 0.2)) +
  # scale_colour_manual(values = rev(met.brewer("Cross", 6))) +
  scale_colour_manual(values = rev(met.brewer("Tiepolo", 6)),
                      breaks = c(100, 250, 500, 1000, 2500, 5000),
                      labels = c("100  ", "250  ", "500  ", "1000  ", "2500  ", "5000  ")) +
  geom_hline(yintercept = c(0.05, 0.8), lty = 2) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9))) +
  theme(legend.text = element_text(size = 18),
        panel.spacing.x = unit(0.75, "cm"),
        panel.spacing.y = unit(0.75, "cm"))

# save png
gg1 +
  theme_png()
ggsave(here("figures", paste(Sys.Date(), "simulations-candidate-gene.png", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height)

# save png without background
gg1 +
  theme_png_no_background()
ggsave(here("figures", paste(Sys.Date(), "simulations-candidate-gene-no-background.png", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height)

# pdf LaTeX
gg1 +
  theme_pdf()
ggsave(here("figures", paste(Sys.Date(), "simulations-candidate-gene.pdf", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height, device = cairo_pdf)


# |-- 2. Tables -----------------------------------------------------------

# variables to filter by
# all possible combinations
dt1_c <- CJ(or = dt1[, unique(or_t)], theta = dt1[, unique(theta0)], n = dt1[, unique(sample_size)])

# get individual datasets and select maximum misclassification when p>=80%
power_p <- function(data, odds, theta, n, p_above) {
  d <- data[or_t == as.numeric(odds) & theta0 == as.numeric(theta) & sample_size == as.numeric(n)]
  values <- ifelse(c(d[p >= p_above, .N] != 0, d[p >= p_above, .N] != 0),
                   c(d[d[p >= p_above, .N], misrate], tail(d[p >= p_above, p], 1)),
                   c(0, max(d[, p])))
  d[1, .(or_t, theta0, sample_size,
         misrate = values[1], p_80 = values[2])]
}

# get data frame
dt1_p80 <- lapply(seq_len(dt1_c[, .N]),
                  function(x) power_p(data = dt1,
                                      odds = as.numeric(dt1_c[x, 1]),
                                      theta   = as.numeric(dt1_c[x, 2]),
                                      n    = as.numeric(dt1_c[x, 3]),
                                      p_above = 0.80)) %>%
  rbindlist()

# save it
fwrite(dt1_p80, here("data", paste(Sys.Date(), "candidate-gene-power-above-80.csv", sep = "_")))

dt1_p80[, or_f := factor(or_t, levels = c(10, 5, 3, 2, 1.5, 1.25))]
dt1_p80[, thetaf := factor(theta0, levels = c(0.5, 0.25, 0.10, 0.05))]

# Table 5
rbind(dcast(dt1_p80[sample_size == 100],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 250],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 500],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 1000], or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 2500], or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 5000], or_f ~ theta0, value.var = "misrate"))


# |--- Just curiosity -----------------------------------------------------

# Pearson's Chi-squared shows differences at high Delta_T and lowest n_i
# According to theory on contingency tables, this makes sense

# how does max (misrate) varies by theta0 across other parameters?
gg_p80 <-
  dt1_p80 %>%
  ggplot(aes(as.factor(theta0), misrate, group = 1)) +
  geom_line() +
  # geom_hline(yintercept = c(0), col = "gray55") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0,1), expand = expansion(c(0,0))) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  # geom_point(size = 2) +
  geom_point(aes(fill = ifelse(misrate != 0, "black", "white")), size = 2, shape = 21) +
  scale_fill_manual(values = c("black", "white")) +
  facet_grid(or_t ~ sample_size,
             labeller = label_bquote(
               rows = Delta[T] == .(as.character(or_t)),
               cols = n[i] == .(as.character(sample_size)))) +
  labs(x = expression(paste("Healthy controls risk factor probability,")~theta[0]),
       y = expression(paste("Misclassification rate,")~gamma))


gg_p80 +
  theme_png() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),
        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted", size = 0.4))
ggsave(here("figures", paste(Sys.Date(), "theta0-inference.png", sep = "_")),
       scale = 0.75, dpi = dpi, width = 14, height = 14)

gg_p80 +
  theme_pdf() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),
        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted", size = 0.4),)
ggsave(here("figures", paste(Sys.Date(), "theta0-inference.pdf", sep = "_")),
       scale = 0.75, dpi = 320, width = 14, height = 14, device = cairo_pdf)



# |- Serology -------------------------------------------------------------

# readjust parameters
width <- 14
height <- 14

# |-- 1. Figures ----------------------------------------------------------

gg2 <-
  dt2 %>%
  ggplot(aes(misrate, p)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6",expression(1-beta), "1")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis",
       colour = expression(paste("Sample size")~(n[0]==n[1])~paste("  "))) +
  facet_grid(sp_f ~ se_f,
             labeller = label_bquote(
               rows = pi[Sp] == .(as.character(sp_f)),
               cols = pi[Se] == .(as.character(se_f)))
  ) +
  geom_line(aes(col = ss_f), lwd = 1.2) +
  # scale_colour_manual(values = darken(viridis(6), 0.2)) +
  # scale_colour_manual(values = rev(met.brewer("Cross", 6))) +
  scale_colour_manual(values = rev(met.brewer("Tiepolo", 6)),
                      breaks = c(100, 250, 500, 1000, 2500, 5000),
                      labels = c("100  ", "250  ", "500  ", "1000  ", "2500  ", "5000  ")) +
  geom_hline(yintercept = c(0.05, 0.8), lty = 2) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9))) +
  theme(legend.text = element_text(size = 18),
        panel.spacing.x = unit(0.75, "cm"),
        panel.spacing.y = unit(0.75, "cm"))

# save png
gg2 +
  theme_png()
ggsave(here("figures", paste(Sys.Date(), "simulations-serology-or-3.png", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height)
# save png without background
gg2 +
  theme_png_no_background()
ggsave(here("figures", paste(Sys.Date(), "simulations-serology-no-background-or-3.png", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height)

# pdf LaTeX
gg2 +
  theme_pdf()
ggsave(here("figures", paste(Sys.Date(), "simulations-serology-or-3.pdf", sep = "_")),
       scale = scale, dpi = dpi, width = width, height = height, device = cairo_pdf)


# |-- 2. Tables -----------------------------------------------------------

# variables to filter by all possible combinations
dt2_c <- CJ(sens = dt2[, unique(sensitivity)], spec = dt2[, unique(specificity)], n = dt2[, unique(sample_size)])

# get individual datasets and select maximum misclassification when p>=80%
power_p <- function(data, sens, spec, n, p_above) {
  d <- data[sensitivity == as.numeric(sens) & specificity == as.numeric(spec) & sample_size == as.numeric(n)]
  values <- ifelse(c(d[p >= p_above, .N] != 0, d[p >= p_above, .N] != 0),
                   c(d[d[p >= p_above, .N], misrate], tail(d[p >= p_above, p], 1)),
                   c(0, max(d[, p])))
  d[1, .(sensitivity, specificity, sample_size,
         misrate = values[1], p_80 = values[2])]
}

# get data frame
dt2_p80 <- lapply(seq_len(dt2_c[, .N]),
                  function(x) power_p(data = dt2,
                                      sens = as.numeric(dt2_c[x, 1]),
                                      spec = as.numeric(dt2_c[x, 2]),
                                      n    = as.numeric(dt2_c[x, 3]),
                                      p_above = 0.80)) %>%
  rbindlist()

# save it
fwrite(dt2_p80, here("data", paste(Sys.Date(), "serology-power-above-80-or-3.csv", sep = "_")))

dt2_p80[, sens := factor(sensitivity, levels = rev(unique(dt2_p80$sensitivity)))]
dt2_p80[, spec := factor(specificity, levels = rev(unique(dt2_p80$specificity)))]

# Table 6
rbind(dcast(dt2_p80[sample_size == 100],  spec ~ sens, value.var = "misrate"),
      dcast(dt2_p80[sample_size == 250],  spec ~ sens, value.var = "misrate"),
      dcast(dt2_p80[sample_size == 500],  spec ~ sens, value.var = "misrate"),
      dcast(dt2_p80[sample_size == 1000], spec ~ sens, value.var = "misrate"),
      dcast(dt2_p80[sample_size == 2500], spec ~ sens, value.var = "misrate"),
      dcast(dt2_p80[sample_size == 5000], spec ~ sens, value.var = "misrate"))









gg2_p80 <-
  dt2_p80 %>%
  ggplot(aes(as.factor(sample_size), misrate, group = 1)) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0,1), expand = expansion(c(0,0))) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  geom_point(aes(fill = ifelse(misrate != 0, "black", "white")), size = 2, shape = 21) +
  scale_fill_manual(values = c("black", "white")) +
  facet_grid(specificity ~ sensitivity,
             labeller = label_bquote(
               rows = pi[sp] == .(as.character(specificity)),
               cols = pi[se] == .(as.character(sensitivity)))) +
  labs(x = expression(paste("Sample size,")~n[i]),
       y = expression(paste("Misclassification rate,")~gamma))


gg2_p80 +
  theme_png() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),
        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted", size = 0.4),)



# Data real-world ---------------------------------------------------------

# data steiner
dt_st <- fread(here("data-cluster/no-correction/2022-05-18_simulation-results-steiner2020.csv"))
dt_st[, snp := factor(snp, levels = c("CTLA4", "PTPN22", "TNF2", "TNF1", "IRF5"))]
dt_st[, signif := ifelse(snp %in% c("CTLA4", "PTPN22"), "yes", "no")]

# data cliff
dt_cl <- fread(here("data-cluster/no-correction/2022-05-18_simulation-results-cliff2019.csv"))
dt_cl[, virus := factor(virus, levels = c("HSV1", "HSV2", "EBV", "CMV", "HHV6", "VZV"))]


# dt_st[p >= 0.8]

cbind(
  dt_st[sensitivity == 1 & specificity == 1 & misrate == 0.24, .(snp, theta0, or_t, theta1, theta1_true)][, 1],
  round(dt_st[sensitivity == 1 & specificity == 1 & misrate == 0.24, .(theta0, or_t, theta1, theta1_true)], 2)
)

cbind(
  dt_cl[sensitivity == 0.975 & specificity == 0.975 & misrate == 0, .(virus)][, 1],
  round(dt_cl[sensitivity == 0.975 & specificity == 0.975 & misrate == 0, .(theta0, or_t, theta1, theta1_true, p)], 2)
)[order(-p)]

# dt_st[misrate == 0.1 & snp == "PTPN22"]
# dt_st[snp == "PTPN22" & p<=0.5]
# dt_cl[virus == "HSV1"]
# dt_cl[virus == "HSV1" & p <= 0.3]
# dt_cl[p > 0.2 & virus != "HSV1"]


# col_line <- met.brewer("Derain")[c(6, 3)]
col_line <- c("#17486f", "#6f9969")
width <- 10
height <- 10


# |- Steiner et al. (2020) ------------------------------------------------

ggst <-
  dt_st %>%
  ggplot(aes(misrate, p)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6",expression(1-beta), "1"),
                     expand = expansion(0,0)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0,1,0.2),
                     labels = c("0","0.2","0.4","0.6","0.8", "1"),
                     expand = expansion(add = c(0.12, 0))) +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Prob. rejecting null hypothesis") +
  geom_line(aes(group = snp, col = signif), lwd = 1.1) +
  geom_text(data = dt_st[misrate == 0],
            aes(label = snp),
            vjust = c(-0.7, 0.2, 1.5, 0.3, 0),
            hjust = c(0.85, 1.05, 1.05, 1.05, 1.05), size = 4) +
  # scale_colour_manual(values = darken(viridis(6), 0.2)) +
  # scale_colour_manual(values = rev(met.brewer("Cross", 6))) +
  scale_colour_manual(values = col_line) +
  geom_hline(yintercept = c(0.05, 0.8), lty = 2) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9))) +
  theme(legend.position = "none",
        panel.spacing.x = unit(0.75, "cm"),
        panel.spacing.y = unit(0.75, "cm"))

ggst +
  theme_png() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-steiner2020.png", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height)
ggst +
  theme_png_no_background() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-steiner2020-no-background.png", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height)
ggst +
  theme_pdf() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-steiner2020.pdf", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height, device = cairo_pdf)


# |-  Cliff et al. (2019) -------------------------------------------------

ggcl <-
  dt_cl %>%
  ggplot(aes(misrate, p)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6",expression(1-beta), "1"),
                     expand = expansion(0,0)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0,1,0.2),
                     labels = c("0","0.2","0.4","0.6","0.8", "1"),
                     expand = expansion(add = c(0.1, 0))) +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis") +
  geom_line(aes(group = virus), lwd = 1.1, colour = col_line[1]) +
  geom_text(data = dt_cl[misrate == 0],
            aes(label = virus),
            vjust = c(0.4,0.2,0.2,0.2,1.1,0.3), hjust = 1.05, size = 4) +
  # scale_colour_manual(values = darken(viridis(6), 0.2)) +
  # scale_colour_manual(values = rev(met.brewer("Cross", 6))) +
  # scale_colour_manual(values = met.brewer("Redon", 6)) +
  geom_hline(yintercept = c(0.05, 0.8), lty = 2) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9))) +
  theme(legend.position = "none",
        panel.spacing.x = unit(0.75, "cm"),
        panel.spacing.y = unit(0.75, "cm"))
ggcl +
  theme_png() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-cliff2019.png", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height)
ggcl +
  theme_png_no_background() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-cliff2019-no-background.png", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height)
ggcl +
  theme_pdf() +
  theme(legend.position = "none")
ggsave(here("figures", paste(Sys.Date(), "simulations-cliff2019.pdf", sep = "_")),
       dpi = dpi, scale = scale, width = width, height = height, device = cairo_pdf)



# |- Figure 3, A and B ----------------------------------------------------

library(patchwork)
# width <- 10
width <- 6

# png
ggst + theme_png() + theme(legend.position = "none") +
  ggcl + labs(y = "") + theme_png() + theme(legend.position = "none") +
  patchwork::plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = c(0.011, 1))
ggsave(here("figures", paste(Sys.Date(), "simulations-real-world.png", sep = "_")), scale = scale, dpi = dpi, width = width*2.2, height = width)


# pdf
ggst + theme_pdf() + theme(legend.position = "none") +
  ggcl + labs(y = "") + theme_pdf() + theme(legend.position = "none") +
  patchwork::plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 25),
        plot.tag.position = c(0.011, 1))
ggsave(here("figures", paste(Sys.Date(), "simulations-real-world.pdf", sep = "_")), scale = 1, dpi = 320, width = width*2.2, height = width, device = cairo_pdf)

# end
