## 56_env_figure.R -- supplementary figure: environment metrics at CO-QTL windows
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
R <- fread(file.path(TAB, "Table44_IBSpy_env_metrics.csv"))
R[, type := factor(type, levels = c("control", "CO-QTL"), labels = c("background", "CO-QTL window"))]
mk <- function(v, lab) {
  p <- wilcox.test(R[[v]] ~ R$type)$p.value
  ggplot(R, aes(x = type, y = get(v), fill = type)) +
    geom_violin(colour = NA, alpha = 0.45, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.22) +
    annotate("text", x = 1.5, y = max(R[[v]], na.rm = TRUE) * 1.02,
             label = sprintf("P = %s", format.pval(p, digits = 2)), size = 2) +
    scale_fill_manual(values = c("#BDBDBD", "#E08A2E"), guide = "none") +
    theme_classic(base_size = 7) + labs(x = NULL, y = lab)
}
p1 <- mk("shannon", "haplotype diversity (Shannon)")
p2 <- mk("mantel_r", "Mantel r (haplotype vs distance)")
p3 <- mk("Fst_lat", "latitudinal differentiation (Fst-like)")
p4 <- mk("lat_slope", "latitude slope of dominant haplotype")
library(patchwork)
p <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "a")
ggsave(file.path(FIG, "Fig4f_IBSpy_env_metrics.png"), p, width = 6.6, height = 4.6, dpi = 600)
cat("env figure done\n")
