suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
img <- function(p) { i <- readPNG(p); g <- rasterGrob(i, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
p2a <- img(file.path(OUT, "panels", "Fig2a_native.png"))
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv")); eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.2, size = 2.1) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.18) +
  theme(legend.position = c(0.78, 0.25), legend.key.size = unit(2.5,"mm"))
p2cd <- img(file.path(OUT, "panels", "Fig2cd_specificity.png"))
p2e <- img(file.path(OUT, "panels", "Fig2c_genes.png"))
fig2 <- p2a / p2b / p2cd / p2e + plot_annotation(tag_levels = c("a","b","c","d")) +
  plot_layout(heights = c(1.0, 0.52, 0.60, 0.80))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 245, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 245, units = "mm", dpi = 600)
cat("Fig2 done\n")
