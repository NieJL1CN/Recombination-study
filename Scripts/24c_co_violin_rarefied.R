## 24c_co_violin_rarefied.R -- Fig3b: per-family CO, raw and rarefied (both kept)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
draw <- function(d, ylab, fname, ymax = NULL) {
  if (!is.null(ymax)) d <- d[co <= ymax]
  ord <- d[, .(m = median(co)), by = family][order(m)]$family
  d[, family := factor(family, levels = ord)]
  p <- ggplot(d, aes(x = family, y = co)) +
    geom_violin(fill = "#2E8B57", colour = NA, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.22) +
    geom_jitter(width = 0.12, size = 0.2, alpha = 0.4, colour = "grey25") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 0.9, colour = "#C0392B") +
    { if (!is.null(ymax)) coord_cartesian(ylim = c(0, ymax)) else NULL } +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5)) +
    labs(x = NULL, y = ylab)
  ggsave(file.path(FIG, fname), p, width = 5.4, height = 3.6, dpi = 600)
}
## raw
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
CO[, co := as.numeric(R1 + R2a + C + R2b + R3)]
draw(CO, "CO per individual (all markers, 725-3,595)", "Fig3b_CO_violin.png")
## rarefied
R <- fread(file.path(TAB, "Table28_CO_rarefied.csv"))
R[, family := sub("_WATDE.*$", "", dataset)]
R[, co := as.numeric(co_rare)]
draw(R, "CO per individual (rarefied to 700 markers)", "Fig3b_CO_violin_rarefied.png", ymax = 90)
## side-by-side summary of the two
A <- CO[, .(family, co, panel = "all markers")]
B <- R[, .(family, co, panel = "700 markers (rarefied)")]
D <- rbind(A, B)
F <- D[, .(mean = mean(co)), by = .(family, panel)]
ord <- F[panel == "all markers"][order(-mean)]$family
F[, family := factor(family, levels = ord)]
p <- ggplot(F, aes(x = family, y = mean, group = family)) +
  geom_segment(aes(xend = family, y = 0, yend = mean), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = panel, shape = panel), size = 1.6) +
  scale_colour_manual(values = c("all markers" = "#C0392B", "700 markers (rarefied)" = "#2C6E9B"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4.4), legend.position = c(0.85, 0.85)) +
  labs(x = NULL, y = "mean CO per family")
ggsave(file.path(FIG, "Fig3b_CO_family_raw_vs_rarefied.png"), p, width = 5.4, height = 3.4, dpi = 600)
cat("raw and rarefied figures written\n")