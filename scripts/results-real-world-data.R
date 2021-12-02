# Joao Malato
# 29/11/2021
# Script for figures on real-world data


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
steiner <- fread(here("data/2021-11-29_simulation-results-steiner2020.csv"))
steiner[, snp := factor(snp, levels = c("CTLA4", "PTPN22", "TNF2", "TNF1", "IRF5"))]

steiner2 <- fread(here("data/2021-11-30_simulation-results-steiner2020-2.csv"))
steiner2[, snp := factor(snp, levels = c("CTLA4", "PTPN22", "TNF2", "TNF1", "IRF5"))]

cliff <- fread(here("data/2021-11-29_simulation-results-cliff2019.csv"))
cliff[, virus := factor(virus)]


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
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "none",
      legend.margin = margin(0, 0, -4, 0),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
)


# Candidate gene ----------------------------------------------------------

# steiner[snp == "CTLA4"] %>%
#   ggplot(aes(misrate, p)) +
#   geom_ribbon(aes(xmax = misrate, xmin = misrate, ymax = p_max, ymin = p_min), alpha = 0.4) +
#   geom_line()

gg_steiner <-
  steiner %>%
  ggplot(aes(misrate, p, group = snp)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  geom_text(data = steiner[misrate == 0],
            aes(x = 0, y = p, label = snp, colour = snp), family = "LM Roman 10",
            hjust = c(0.75, 0.9, 1.1, 1.1, 1.1),
            vjust = c(-0.4, -0.4, 0.4, -0.4, -0.6), alpha = 1, size = 4.5, show.legend = FALSE) +
  geom_line(aes(colour = snp), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(5), 0.2)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1"), expand = expansion(0,0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1"), expand = expansion(add = c(0.07, 0))) +
  coord_cartesian(expand = TRUE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis")
gg_steiner
# png
width <- 10
ggsave(here("figures/simulations-steiner2020.png"), gg_steiner, scale = 1, dpi = 320, width = width, height = width*0.6)
# pdf LaTeX
ggsave(here("figures/simulations-steiner2020.pdf"), gg_steiner, scale = 1, dpi = 320, width = width, height = width*0.6, device = cairo_pdf)


gg_steiner2 <-
  steiner2 %>%
  ggplot(aes(misrate, p, group = snp)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  geom_text(data = steiner2[misrate == 0],
            aes(x = 0, y = p, label = snp, colour = snp), family = "LM Roman 10",
            hjust = 1.1, vjust = c(0.1, 0.1, 0.65, -0.1, -0.3), alpha = 1, size = 4.5, show.legend = FALSE) +
  geom_line(aes(colour = snp), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(5), 0.2)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1"), expand = expansion(0,0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1"), expand = expansion(add = c(0.07, 0))) +
  coord_cartesian(expand = TRUE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis")
gg_steiner2


# Serology study ----------------------------------------------------------

gg_cliff <-
  cliff %>%
  ggplot(aes(misrate, p, group = virus)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  geom_text(data = cliff[misrate == 0],
            aes(x = 0, y = p, label = virus, colour = virus), family = "LM Roman 10",
            hjust = 1.1, vjust = c(0.2, 0.2, 0.3, 0.2, 1.5, 0.1), alpha = 1, size = 4.5, show.legend = FALSE) +
  geom_line(aes(colour = virus), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(6), 0.2)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1"), expand = expansion(0,0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1"), expand = expansion(add = c(0.07, 0))) +
  coord_cartesian(expand = TRUE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis")
gg_cliff
width <- 10
# png
ggsave(here("figures/simulations-cliff2019.png"), gg_cliff, scale = 1, dpi = 320, width = width, height = width*0.6)
# pdf LaTeX
ggsave(here("figures/simulations-cliff2019.pdf"), gg_cliff, scale = 1, dpi = 320, width = width, height = width*0.6, device = cairo_pdf)



# Both plots --------------------------------------------------------------

library(patchwork)

gg_steiner / gg_cliff +
  patchwork::plot_annotation(tag_levels = "A")

#
# dt_m <- rbind(steiner[, .(misrate, p, var = snp, plot = "steiner")], cliff[, .(misrate, p, var = virus, plot = "cliff")])
#
# dt_m %>%
#   ggplot(aes(misrate, p)) +
#   facet_wrap(~plot) +
#   annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
#   geom_text(data = dt_m[misrate == 0],
#             aes(x = 0, y = p, label = var, colour = var), family = "LM Roman 10",
#             hjust = 1.1, vjust = 0.1, alpha = 1, size = 4.5, show.legend = FALSE) +
#   geom_line(aes(colour = var), size = 1.1) +
#   geom_hline(yintercept = 0.05, alpha = 0.8) +
#   # scale_colour_manual(values = darken(viridis(12), 0.2)) +
#   scale_y_continuous(limits = c(0,1),
#                      breaks = sort(c(0.05, seq(0,1, 0.2))),
#                      labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1"), expand = expansion(0,0)) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1"), expand = expansion(add = c(0.07, 0))) +
#   coord_cartesian(expand = TRUE, clip = "off") +
#   labs(x = expression(paste("Misclassification rate,")~gamma),
#        y = "Probability of rejecting null hypothesis")
