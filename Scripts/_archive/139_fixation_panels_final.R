## 139_fixation_panels_final.R -- single panels, one-row legend, both metrics
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
draw <- function(dat, yvar, ylab, tag) {
  CR <- dat[order(gen), .(cross = { i <- which(.SD[[1]] >= 0.5); if (length(i)) gen[i[1]] else NA_integer_ }),
            by = .(lab), .SDcols = yvar]
  CR <- CR[!is.na(cross)]
  p <- ggplot(dat, aes(gen, .data[[yvar]], colour = lab)) +
    geom_line(aes(linewidth = lab == "No QTL")) +
    geom_point(data = CR, aes(cross, 0.5), size = 0.55, show.legend = FALSE) +
    geom_text(data = CR, aes(cross, 0.42, label = cross), size = 1.5, show.legend = FALSE) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
    scale_linewidth_manual(values = c("TRUE" = 0.45, "FALSE" = 0.2), guide = "none") +
    scale_colour_manual(values = pal, name = NULL) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_continuous(breaks = 1:5, limits = c(1, 5)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
    labs(x = "Generation", y = ylab) +
    theme_classic(base_size = 6.5) +
    theme(legend.position = "bottom", legend.key.width = unit(7, "pt"),
          legend.text = element_text(size = 4.6), legend.spacing.x = unit(1.5, "pt"),
          legend.margin = margin(0, 0, 0, 0), axis.text = element_text(size = 5.5),
          axis.title = element_text(size = 6))
  ggsave(file.path(ROOT, sprintf("Figures/FigS_%s_%dcM.png", tag, d)), p,
         width = 4.4, height = 2.3, dpi = 600)
}
for (d in c(1, 2, 5, 10)) {
  draw(S[dist == d], "TTdd", "TTdd (target homozygous, drag removed)", "fixation")
  draw(S[dist == d], "cleanT", "T carriers free of the deleterious allele", "breakage")
}
cat("done\n")