## 149_per_family_figures.R -- one fixation and one breakage figure per family
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table127_sim_tight_distances.csv"))
fam <- data.table(family = c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
                  pct = c(17.3, 13.9, 12.3, 6.6, 4.5, 3.5, 2.3),
                  label = c("ParW049 (7B, +17.3%)", "ParW579 (3B+4B, +13.9%)", "ParW729 (1B+4A, +12.3%)",
                            "RobW145 (6B, +6.6%)", "ParW007 (6A, +4.5%)", "ParW352 (7D, +3.5%)",
                            "ParW145 (5D, +2.3%)"))
pal <- c("0.1 cM" = "#7B241C", "0.5 cM" = "#C0392B", "1 cM" = "#E08A2E",
         "2 cM" = "#2E8B57", "3 cM" = "#2C6E9B")
draw <- function(dat, famlab, yvar, ylab, tag) {
  p <- ggplot(dat, aes(gen, .data[[yvar]], colour = dist_f, linetype = type,
                       group = interaction(dist_f, type))) +
    geom_line(linewidth = 0.42) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
    scale_colour_manual(values = pal, name = "Target-drag distance") +
    scale_linetype_manual(values = c("No QTL" = "solid", "With QTL" = "22"),
                          labels = c("No QTL", famlab), name = NULL) +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
    labs(x = "Generation", y = ylab) +
    theme_classic(base_size = 7) +
    theme(legend.position = "right", legend.key.size = unit(8, "pt"),
          legend.text = element_text(size = 5.2), legend.key.width = unit(14, "pt"))
  ggsave(file.path(ROOT, sprintf("Figures/FigS_%s_%s.png", tag, famlab_short)), p,
         width = 4.0, height = 2.5, dpi = 600)
}
for (i in seq_len(nrow(fam))) {
  famlab_short <- fam$family[i]
  famlab <- fam$label[i]
  D <- R[family %in% c("baseline", fam$family[i]),
         .(TTdd = mean(TTdd, na.rm = TRUE), cleanT = mean(cleanT, na.rm = TRUE)),
         by = .(family, dist, gen)]
  D[, type := fifelse(family == "baseline", "No QTL", "With QTL")]
  D[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(0.1,0.5,1,2,3), " cM"))]
  draw(D, famlab, "TTdd", "TTdd (target homozygous, drag removed)", "fixation")
  draw(D, famlab, "cleanT", "T carriers free of the deleterious allele", "breakage")
  cat("done", fam$family[i], "\n")
}
