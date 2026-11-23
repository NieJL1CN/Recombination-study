suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
q <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                chr = c("3B","1B","4A","7D","4B","7B","6A","5D","6B"),
                Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8),
                trait = c("distal","distal","distal","distal","distal","centromere","whole","whole","whole"))
len[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
q[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
len[, cen := Mb * 0.5]
p2a <- ggplot() +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = 0, ymax = Mb), fill = "grey85", colour = "grey35", linewidth = 0.25) +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = cen - 4, ymax = cen + 4), fill = "white", colour = NA) +
  geom_point(data = q, aes(idx, Mb, colour = trait), size = 1.1) +
  geom_text_repel(data = q, aes(idx, Mb, label = qtl), size = 2.0, direction = "y", nudge_x = 0.45,
                  segment.size = 0.15, min.segment.length = 0) +
  scale_colour_manual(values = c(distal = "#2E6FA3", centromere = "grey30", whole = "#E08214")) +
  scale_x_continuous(breaks = seq_along(chrlev), labels = chrlev, expand = expansion(add = 0.6)) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = NULL, y = "Physical position (Mb)", colour = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5,"mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
ggsave(file.path(OUT, "panels", "Fig2a_native.png"), p2a, width = 150, height = 85, units = "mm", dpi = 600)
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.2, size = 2.3) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.18)
p2a_img <- { img <- readPNG(file.path(OUT, "panels", "Fig2a_native.png")); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
fig2 <- p2a_img / p2b + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.15, 1))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 185, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 185, units = "mm", dpi = 600)
cat("Fig2 a+b done\n")
