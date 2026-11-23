suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(cowplot)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"), axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9, face = "bold"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7), plot.margin = margin(2, 3, 2, 2)))
img_panel <- function(path) {
  img <- readPNG(path); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void()
}
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
mb <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                 Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8))
eff <- merge(eff, mb, by = "qtl")
p2a <- img_panel(file.path(OUT, "panels", "metaQTL_ideogram.png"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 0.6) + geom_point(aes(y = mean_A), size = 0.4) + geom_point(aes(y = mean_B), size = 0.4) +
  coord_flip() + scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL) + theme(legend.position = "none")
fig2 <- (p2a | p2b) / p2c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 170, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 170, units = "mm", dpi = 600)
cat("Fig2 rebuilt\n")
