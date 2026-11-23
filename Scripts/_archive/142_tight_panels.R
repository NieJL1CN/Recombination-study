## 142_tight_panels.R -- panels for 0.1, 0.5, 1, 2, 3 cM
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table127_sim_tight_distances.csv"))
S <- R[, .(TTdd = mean(TTdd, na.rm = TRUE)), by = .(family, pct, dist, gen)]
labs <- c(baseline = "No QTL", ParW049 = "ParW049", ParW579 = "ParW579 (3B+4B)",
          ParW729 = "ParW729 (1B+4A)", RobW145 = "RobW145", ParW007 = "ParW007",
          ParW352 = "ParW352", ParW145 = "ParW145")
S[, lab := labs[family]]
S[, lab := factor(lab, levels = c("No QTL","ParW049","ParW579 (3B+4B)","ParW729 (1B+4A)",
                                  "RobW145","ParW007","ParW352","ParW145"))]
pal <- c("No QTL" = "grey20", "ParW049" = "#C0392B", "ParW579 (3B+4B)" = "#2C6E9B",
         "ParW729 (1B+4A)" = "#2E8B57", "RobW145" = "#6A4C93", "ParW007" = "#E08A2E",
         "ParW352" = "#8C564B", "ParW145" = "#7F7F7F")
draw <- function(dat, tag, d, xmax) {
  p <- ggplot(dat, aes(gen, TTdd, colour = lab)) +
    geom_line(aes(linewidth = lab == "No QTL")) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
    scale_linewidth_manual(values = c("TRUE" = 0.45, "FALSE" = 0.2), guide = "none") +
    scale_colour_manual(values = pal, name = NULL) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_continuous(breaks = 1:xmax, limits = c(1, xmax)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
    labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
    theme_classic(base_size = 6.5) +
    theme(legend.position = "bottom", legend.key.width = unit(7, "pt"),
          legend.text = element_text(size = 4.6), legend.spacing.x = unit(1.5, "pt"),
          axis.text = element_text(size = 5.5), axis.title = element_text(size = 6))
  ggsave(file.path(ROOT, sprintf("Figures/FigS_%s_%scM.png", tag, d)), p,
         width = 4.4, height = 2.3, dpi = 600)
}
for (d in c(0.1, 0.5, 1, 2, 3)) {
  draw(S[dist == d], "fixation", d, if (d >= 2) 5 else 10)
}
## baseline comparison
B <- S[family == "baseline", .(TTdd = TTdd), by = .(dist, gen)]
B[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(0.1,0.5,1,2,3), " cM"))]
p <- ggplot(B, aes(gen, TTdd, colour = dist_f)) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.55) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = c("0.1 cM" = "#7B241C", "0.5 cM" = "#C0392B", "1 cM" = "#E08A2E",
                                 "2 cM" = "#2E8B57", "3 cM" = "#2C6E9B"), name = "Target-drag distance") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_baseline_distances.png"), p, width = 3.8, height = 2.4, dpi = 600)
cat("done\n")
