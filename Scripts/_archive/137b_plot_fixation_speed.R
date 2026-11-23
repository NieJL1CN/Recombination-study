## 137b_plot_fixation_speed.R -- fixation and linkage-breakage curves
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table123_sim_trajectories.csv"))
S <- R[, .(TTdd = mean(TTdd, na.rm = TRUE), cleanT = mean(cleanT, na.rm = TRUE)),
       by = .(family, pct, dist, gen)]
S[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(1,2,5,10), " cM"))]
S[, fam_lab := fifelse(family == "baseline", "No QTL", paste0(family, " (+", pct, "%)"))]
lv <- c("No QTL", paste0(c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
                        " (+", c(17.3,13.9,12.3,6.6,4.5,3.5,2.3), "%)"))
S[, fam_lab := factor(fam_lab, levels = lv)]
pal <- c("No QTL" = "grey25", "ParW049 (+17.3%)" = "#C0392B", "ParW579 (+13.9%)" = "#2C6E9B",
         "ParW729 (+12.3%)" = "#2E8B57", "RobW145 (+6.6%)" = "#6A4C93",
         "ParW007 (+4.5%)" = "#E08A2E", "ParW352 (+3.5%)" = "#8C564B",
         "ParW145 (+2.3%)" = "#7F7F7F")
mk <- function(y, ylab, ttl) {
  ggplot(S, aes(gen, .data[[y]], colour = fam_lab)) +
    geom_line(aes(linewidth = fam_lab == "No QTL", alpha = fam_lab == "No QTL")) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.22, colour = "grey45") +
    facet_wrap(~ dist_f, nrow = 1) +
    scale_linewidth_manual(values = c("TRUE" = 0.8, "FALSE" = 0.35), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.7), guide = "none") +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
    labs(x = if (y == "TTdd") NULL else "Generation of selfing with marker-assisted selection",
         y = ylab, title = ttl) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6.5),
          plot.title = element_text(size = 6.8), legend.position = "none")
}
p1 <- mk("TTdd", "TTdd (target homozygous, drag removed)", "a  Fixation of the clean target haplotype")
p2 <- mk("cleanT", "T carriers free of the deleterious allele", "b  Breakage of the target-drag linkage")
p <- p1 / p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.key.size = unit(6.5, "pt"),
        legend.text = element_text(size = 5.2), legend.box = "horizontal")
ggsave(file.path(ROOT, "Figures/FigS_fixation_speed.png"), p, width = 7.2, height = 4.4, dpi = 600)
cat("figure done\n")
