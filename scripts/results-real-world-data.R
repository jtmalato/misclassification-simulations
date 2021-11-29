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
steiner[, snp := factor(snp, levels = rev(c("CTLA4", "PTPN22", "TNF2", "TNF1", "IRF5")))]

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
      legend.background = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.spacing.x = unit(3, 'pt'),
      legend.margin = margin(0, 0, -4, 0),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
)


# Candidate gene ----------------------------------------------------------

gg_steiner <-
  steiner %>%
  ggplot(aes(misrate, p, group = snp)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.80, ymax=1, alpha=0.2, fill="gray35") +
  geom_text(data = steiner[misrate == 0 & snp != "IRF5"],
            aes(x = 0, y = p, label = snp, colour = snp),
            hjust = 1.1, vjust = 0.1, alpha = 1, size = 4, show.legend = FALSE) +
  geom_text(data = steiner[misrate == 0 & snp == "IRF5"],
            aes(x = 0, y = p, label = snp, colour = snp),
            hjust = 1.1, vjust = 0.8, alpha = 1, size = 4, show.legend = FALSE) +
  geom_line(aes(colour = snp), size = 1.1) +
  geom_hline(yintercept = 0.05, alpha = 0.8) +
  scale_colour_manual(values = darken(viridis(5), 0.2), guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(limits = c(0,1),
                     breaks = sort(c(0.05, seq(0,1, 0.2))),
                     labels = c("", expression(alpha), "0.2","0.4","0.6","0.8", "1"), expand = expansion(0,0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2), labels = c("0","0.2","0.4","0.6","0.8", "1"), expand = expansion(add = c(0.07, 0))) +
  coord_cartesian(expand = TRUE, clip = "off") +
  labs(x = expression(paste("Misclassification rate,")~gamma),
       y = "Probability of rejecting null hypothesis",
       colour = "SNP") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size=2, alpha = 0.9)))
gg_steiner
# png
ggsave(here("figures/simulations-steiner2020.png"), gg_steiner,
       scale = 1, dpi = 320, width = 13, height = 9)
# pdf LaTeX
ggsave(here("figures/simulations-steiner2020.pdf"), gg_steiner,
       scale = 1, dpi = 320, width = 13, height = 9, device = cairo_pdf)

