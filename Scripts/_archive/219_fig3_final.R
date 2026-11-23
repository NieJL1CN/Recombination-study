suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(png); library(grid)})
OUT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Figures/main"
img <- function(p) { i <- readPNG(p); g <- rasterGrob(i, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
a1 <- img(file.path(OUT, "panels", "Fig3a_coloc.png"))
a2 <- img(file.path(OUT, "panels", "Fig3a2_metaqlt.png"))
b <- img(file.path(OUT, "panels", "Fig3b_profiles.png"))
c <- img(file.path(OUT, "panels", "Fig3c_within_family.png"))
fig <- (a1 | a2) / (b | c) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.15, 0.75), widths = c(1.25, 0.85))
ggsave(file.path(OUT, "Fig3.tiff"), fig, width = 183, height = 180, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig, width = 183, height = 180, units = "mm", dpi = 600)
cat("Fig3 composed\n")
