## ============================================================================
## 02_map_figures.R -- marker density (CMplot + ideogram), linkage-map view,
## collinearity (3 x 7 facets)
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER  <- file.path(ROOT, "Source/Derived"); FIG <- file.path(ROOT, "Figures")
CHR  <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COL  <- c(blue = "#2C6E9B", purple = "#6A4C93", green = "#2E8B57", orange = "#E08A2E")
PAL  <- c("AXIOM.1" = COL[["blue"]], "AXIOM.2" = COL[["purple"]], "KASP" = COL[["green"]])
L <- readRDS(file.path(DER, "maps.rds")); M <- L$markers; CM <- L$consensus
pos <- unique(M[, .(marker, ppos)])
setwd(FIG)

for (P in unique(CM$panel)) {
  D <- merge(CM[panel == P, .(chr, marker, cM)], pos, by = "marker")
  if (nrow(D) < 100) next

  ## CMplot density is produced by 02b_cmplot.R (bp units, Mb labels)

  ## ---- 1b. ideogram with green marker ticks -----------------------------
  D[, chr := factor(chr, levels = rev(CHR))]
  ends <- D[, .(xend = max(ppos) / 1e6), by = chr]
  p <- ggplot(D, aes(x = ppos / 1e6, y = chr)) +
    geom_segment(data = ends, aes(x = 0, xend = xend, y = chr, yend = chr),
                 colour = "grey90", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(size = 0.16, colour = COL[["green"]], alpha = 0.8) +
    theme_classic(base_size = 7) +
    theme(axis.text.y = element_text(size = 5, colour = "grey20"),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_blank()) +
    labs(x = "physical position (Mb)", y = NULL)
#   ggsave(file.path(FIG, paste0("Fig1b_ideogram_", P, ".png")), p, width = 5, height = 3.4, dpi = 500)   ## removed by request

  ## ---- 1b. LinkageMapView-style genetic map -----------------------------
  G <- D
  G[, chr := factor(chr, levels = rev(CHR))]
  G[, len := max(cM), by = chr]
  G[, y := as.numeric(chr)]
  tick <- G[, .(pos = seq(0, max(cM), by = ifelse(max(cM) > 200, 50, 20))), by = chr]
  tick[, y := match(chr, levels(G$chr))]
  p2 <- ggplot(G, aes(x = cM, y = y)) +
    geom_segment(aes(x = 0, xend = len, yend = y), colour = "grey25", linewidth = 0.45) +
    geom_segment(aes(x = cM, xend = cM, yend = y - 0.18, y = y - 0.18), linewidth = 0.2, colour = COL[["blue"]]) +
    geom_text(data = tick, aes(x = pos, y = y - 0.3, label = pos), size = 1.6, colour = "grey35") +
    scale_y_continuous(breaks = seq_along(levels(G$chr)), labels = levels(G$chr)) +
    theme_classic(base_size = 7) +
    theme(axis.text.y = element_text(size = 5, colour = "grey20"),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_blank()) +
    labs(x = "genetic position (cM)", y = NULL)
#   ggsave(file.path(FIG, paste0("Fig1b_linkageMapView_", P, ".png")), p2, width = 5.6, height = 4, dpi = 500)   ## removed by request

  ## ---- 1d. collinearity 3 x 7 -------------------------------------------
  F <- copy(D)
  F[, grp := factor(substr(chr, 2, 2), levels = c("A", "B", "D"))]
  F[, idx := as.integer(substr(chr, 1, 1))]
  p3 <- ggplot(F, aes(x = ppos / 1e6, y = cM)) +
    geom_point(size = 1.1, alpha = 0.7, colour = PAL[P], stroke = 0) +
    facet_grid(grp ~ idx, scales = "free", switch = "y") +
    theme_bw(base_size = 5.5) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = 5.5, colour = "grey20"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.15, colour = "grey95"),
          panel.border = element_rect(linewidth = 0.25, colour = "grey45"),
          axis.ticks = element_line(linewidth = 0.2), axis.text = element_text(size = 4.6),
          plot.margin = margin(2, 2, 2, 2)) +
    labs(x = "physical position (Mb)", y = "genetic position (cM)")
  ggsave(file.path(FIG, paste0("Fig1d_collinearity_", P, ".png")), p3, width = 4.6, height = 2.3, dpi = 600)
}
cat("map figures done\n")