## 33_co_qtl_figures.R -- Manhattan plots for the four CO phenotypes (own axes per panel)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
P <- fread(file.path(TAB, paste0("Table20_CO_QTL_profiles", Sys.getenv("SUFFIX", ""), ".csv")))
Q <- fread(file.path(TAB, paste0("Table20b_CO_QTL_peaks", Sys.getenv("SUFFIX", ""), ".csv")))
P[, chr := factor(chr, levels = CHR)]
P[, cum := as.numeric(chr) - 0.5 + (pos - min(pos)) / max(1, diff(range(pos))) * 0.9, by = .(dataset, trait, chr)]
fams <- unique(P$dataset)
trait_lab <- c(distal = "distal (R1+R3)", middle = "middle (R2a+R2b)",
               centromere = "centromere (C)", whole = "whole genome")
for (d in fams) {
  D <- P[dataset == d]; if (!nrow(D)) next
  D[, sig := lod > thr]
  ymax <- max(c(D$lod, D$thr), na.rm = TRUE) * 1.05
  traits <- intersect(c("distal", "middle", "centromere", "whole"), unique(D$trait))
  plist <- list()
  for (tr in traits) {
    Di <- D[trait == tr][order(cum)]
    if (!nrow(Di)) { plist[[length(plist) + 1]] <- plot_spacer(); next }
    plist[[length(plist) + 1]] <- ggplot(Di, aes(x = cum, y = lod)) +
      geom_line(colour = "grey70", linewidth = 0.15) +
      geom_point(data = Di[sig == FALSE], colour = "grey60", size = 0.25, stroke = 0) +
      geom_point(data = Di[sig == TRUE], colour = "#C0392B", size = 0.4, stroke = 0) +
      geom_hline(yintercept = unique(Di$thr)[1], colour = "grey35", linetype = "dashed", linewidth = 0.2) +
      scale_x_continuous(breaks = seq_along(CHR), labels = CHR, expand = c(0.01, 0.01)) +
      coord_cartesian(ylim = c(0, ymax)) +
      theme_classic(base_size = 5) +
      theme(axis.text.x = element_text(size = 3, angle = 90, vjust = 0.5),
            axis.text.y = element_text(size = 3.2),
            axis.line = element_line(linewidth = 0.22),
            axis.ticks = element_line(linewidth = 0.18),
            plot.title = element_text(size = 4.2, hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1)) +
      labs(x = NULL, y = "LOD", title = trait_lab[[tr]])
  }
  p <- wrap_plots(plist, ncol = length(traits), nrow = 1)
  ggsave(file.path(FIG, paste0("Fig3g_CO_QTL", Sys.getenv("SUFFIX", ""), "_", sub("_WATDE.*$", "", d), ".png")),
         p, width = 1.5 * length(traits) + 0.4, height = 1.8, dpi = 600)
}
cat("figures for", length(fams), "families\n")
