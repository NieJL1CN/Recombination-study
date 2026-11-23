## 138_fixation_panels.R -- single-panel fixation figures (no compositing)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table123_sim_trajectories.csv"))
S <- R[, .(TTdd = mean(TTdd, na.rm = TRUE), cleanT = mean(cleanT, na.rm = TRUE)),
       by = .(family, pct, dist, gen)]
S <- S[gen <= 5]
labs <- c(baseline = "No QTL", ParW049 = "ParW049", ParW579 = "ParW579 (3B+4B)",
          ParW729 = "ParW729 (1B+4A)", RobW145 = "RobW145", ParW007 = "ParW007",
          ParW352 = "ParW352", ParW145 = "ParW145")
S[, lab := labs[family]]
S[, lab := factor(lab, levels = c("No QTL","ParW049","ParW579 (3B+4B)","ParW729 (1B+4A)",
                                  "RobW145","ParW007","ParW352","ParW145"))]
pal <- c("No QTL" = "grey20", "ParW049" = "#C0392B", "ParW579 (3B+4B)" = "#2C6E9B",
         "ParW729 (1B+4A)" = "#2E8B57", "RobW145" = "#6A4C93", "ParW007" = "#E08A2E",
         "ParW352" = "#8C564B", "ParW145" = "#7F7F7F")
for (d in c(1, 2, 5, 10)) {
  D <- S[dist == d]
  CR <- D[order(gen), .(cross = { i <- which(TTdd >= 0.5); if (length(i)) gen[i[1]] else NA_integer_ }),
          by = .(lab)]
  CR <- CR[!is.na(cross)]
  p <- ggplot(D, aes(gen, TTdd, colour = lab)) +
    geom_line(aes(linewidth = lab == "No QTL")) +
    geom_point(data = CR, aes(cross, 0.5), size = 0.6, show.legend = FALSE) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.2, colour = "grey50") +
    scale_linewidth_manual(values = c("TRUE" = 0.55, "FALSE" = 0.22), guide = "none") +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = 1:5, limits = c(1, 5)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
    labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
    theme_classic(base_size = 7) +
    theme(legend.position = "bottom", legend.key.width = unit(9, "pt"),
          legend.text = element_text(size = 5.2), legend.box = "horizontal",
          legend.spacing.x = unit(2, "pt"))
  ggsave(file.path(ROOT, sprintf("Figures/FigS_fixation_%dcM.png", d)), p,
         width = 3.0, height = 2.6, dpi = 600)
}
cat("panels done\n")
