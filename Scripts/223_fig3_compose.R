suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(png); library(grid)})
OUT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Figures/main"
img <- function(p) { i <- readPNG(p); g <- rasterGrob(i, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
a <- img(file.path(OUT, "panels", "Fig3a_integrated.png"))
b <- img(file.path(OUT, "panels", "Fig3b_correlation.png"))
c <- img(file.path(OUT, "panels", "Fig3b_profiles.png"))
fig <- a / (b | c) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.5, 0.7), widths = c(0.9, 1.3))
ggsave(file.path(OUT, "Fig3.tiff"), fig, width = 183, height = 185, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig, width = 183, height = 185, units = "mm", dpi = 600)
cat("Fig3 composed\n")
