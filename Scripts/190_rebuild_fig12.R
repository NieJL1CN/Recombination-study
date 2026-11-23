suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
img_panel <- function(path) { img <- readPNG(path); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
## Fig1
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1a <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.7,0.8))
p1b <- img_panel(file.path(OUT, "panels", "Fig1b_native.png"))
regm <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1c <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.2) + labs(x = NULL, y = "Crossovers per RIL")
chr <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
F <- merge(chr, len, by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 500, y = max(F$co)*0.95, label = "r = 0.76", size = 2.4, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
design <- "ab
cb
de"
fig1 <- p1a + p1b + p1c + p1d + p1e + plot_layout(design = design, widths = c(2, 1), heights = c(1,1,1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 190, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 190, units = "mm", dpi = 600)
## Fig2
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2a <- img_panel(file.path(OUT, "panels", "Fig2a_native.png"))
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 0.6) + geom_point(aes(y = mean_A), size = 0.4) + geom_point(aes(y = mean_B), size = 0.4) +
  coord_flip() + scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL) + theme(legend.position = "none")
fig2 <- p2a / (p2b | p2c) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.2, 1))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 190, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 190, units = "mm", dpi = 600)
cat("done\n")
