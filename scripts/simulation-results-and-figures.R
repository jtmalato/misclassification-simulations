# Joao Malato
# 26/11/2021
# Script for figures and tables


# Libraries ---------------------------------------------------------------
extrafont::loadfonts(device="win")
library(data.table)
library(here)
library(magrittr)
library(ggplot2)
library(lemon)
library(viridis)
library(colorspace)

# Data --------------------------------------------------------------------
dt1 <- fread(here("data/2021-11-26_simulation-results-candidate-gene.csv"))
dt1[, ss := factor(sample_size)]
dt1[, or_f := factor(or_t, levels = rev(as.character(unique(or_t))))]


# Theme -------------------------------------------------------------------
theme_set(
    theme_classic(base_size = 18, base_family = "LM Roman 10") +
    theme(
      panel.border = element_blank(),
      # plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted", size = 0.4),
      panel.grid.major.y = element_line(colour = c("gray", NA, rep("gray", 5)), linetype = "dotted", size = 0.4),
      axis.title.x = element_text(vjust = -3.5),
      axis.title.y = element_text(vjust = 4),
      strip.background = element_blank(),
      strip.text = element_text(size = 20),
      plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"),
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
)


# Candidate gene ----------------------------------------------------------

# Figure ------------------------------------------------------------------
dt1 %>%
  ggplot(aes(misrate, p, group = ss)) +
  annotate("rect", xmin=0, xmax=1, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  facet_grid(or_f ~ theta0, labeller = label_bquote(rows = Delta[T] == .(as.character(or_f)), cols = theta[0] == .(theta0))) +
  geom_line(aes(colour = ss), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(6), 0.1)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(seq(0,1, 0.2))),
                     labels = c("", "0.2","0.4","0.6","0.8", "1"), minor_breaks = 0.05)

gg1 <- dt1 %>%
  ggplot(aes(misrate, p, group = ss)) +
  annotate("rect", xmin=0, xmax=1, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  facet_grid(or_f ~ theta0, labeller = label_bquote(rows = Delta[T] == .(as.character(or_f)), cols = theta[0] == .(theta0))) +
  geom_line(aes(colour = ss), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(6), 0.1)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis",
       colour = expression(paste("Sample size")~(n[0]==n[1]))) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9)))
# png
ggsave(here("figures/simulations-candidate-gene.png"), gg1,
       scale = 1, dpi = 320, width = 13, height = 14)
# pdf LaTeX
ggsave(here("figures/simulations-candidate-gene.pdf"), gg1,
       scale = 1, dpi = 320, width = 13, height = 14, device = cairo_pdf)


dt1 %>%
  ggplot(aes(misrate, overall_or)) +
  facet_grid(or_f ~ theta0, labeller = label_bquote(rows = Delta[T] == .(as.character(or_f)), cols = theta[0] == .(theta0)), scales = "free_y") +
  geom_line()


# Table -------------------------------------------------------------------

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
# setorder(table_p_80, n, theta0, -or_t, misrate)
fwrite(dt1_p80, here("data/candidate-gene-power-above-80.csv"))

dt1_p80[, or_f := factor(or_t, levels = c(10, 5, 2, 1.5, 1.25))]
dt1_p80[, thetaf := factor(theta0, levels = c(0.5, 0.25, 0.10, 0.05))]

rbind(dcast(dt1_p80[sample_size == 100],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 250],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 250],  or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 1000], or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 2500], or_f ~ theta0, value.var = "misrate"),
      dcast(dt1_p80[sample_size == 5000], or_f ~ theta0, value.var = "misrate"))

# Serology test -----------------------------------------------------------


# Figure ------------------------------------------------------------------

# Table -------------------------------------------------------------------


# end
