## ============================================================================
## 07_qtl_figures.R -- Manhattan plots, QTL summary, meta "snake" map
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
PAL <- c("consensus" = "#2C6E9B", "original" = "#6A4C93")
TRC <- c("PH_M_cm" = "#2C6E9B", "GY_Calc_tha" = "#2E8B57", "GW_M_g1000grn" = "#E08A2E")
PK <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))
PR <- fread(file.path(TAB, "Table2d_QTL_profiles.csv"))

## ---- 1. summary table ------------------------------------------------------
S <- PK[, .(n_qtl = .N, max_lod = round(max(lod), 2), mean_r2 = round(mean(r2), 3)),
        by = .(dataset, trait, map, panel)]
CM <- dcast(S, dataset + trait + panel ~ map, value.var = c("n_qtl", "max_lod"))
fwrite(CM, file.path(TAB, "Table2e_QTL_consensus_vs_original_summary.csv"))
cat("families with QTL:", uniqueN(PK$dataset), "\n")
print(PK[, .(n = .N, med_lod = round(median(lod), 2), max_lod = round(max(lod), 2)), by = map])

## ---- 2. Manhattan per family ----------------------------------------------
PR[, chr := factor(chr, levels = CHR)]
PR[, cum := {
  rg <- diff(range(pos))
  frac <- if (!is.finite(rg) || rg <= 0) rank(pos) / .N else (pos - min(pos)) / rg
  as.numeric(chr) - 0.5 + frac * 0.9
}, by = .(dataset, trait, map, chr)]
fams <- unique(PR$dataset)
suppressPackageStartupMessages(library(patchwork))
for (d in fams) {
  D <- PR[dataset == d]
  if (!nrow(D)) next
  D[, sig := lod > thr]
  ymax <- max(c(D$lod, D$thr), na.rm = TRUE) * 1.05
  maps <- unique(D$map); traits <- unique(D$trait)
  plist <- list()
  for (mp in maps) for (tr in traits) {
    Di <- D[map == mp & trait == tr]
    if (!nrow(Di)) { plist[[length(plist) + 1]] <- plot_spacer(); next }
    plist[[length(plist) + 1]] <- ggplot(Di[order(cum)], aes(x = cum, y = lod)) +
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
      labs(x = NULL, y = "LOD", title = paste0(tr, "  (", mp, ")"))
  }
  p <- wrap_plots(plist, ncol = length(traits), nrow = length(maps), byrow = TRUE)
  ggsave(file.path(FIG, paste0("Fig2a_manhattan_", sub("_WATDE.*$", "", d), ".png")),
         p, width = 1.5 * length(traits) + 0.4, height = 1.35 * length(maps) + 0.3, dpi = 600)
}
# snake map removed: superseded by 11_meta_map.R