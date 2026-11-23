suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
p2a_img <- { img <- readPNG(file.path(OUT, "panels", "Fig2a_native.png")); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.2, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.18)
genes <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
cand <- fread(file.path(TAB, "Table58_candidate_gene_report.csv"))
ord <- unique(genes[order(chr, peak_Mb), qtl_id])
genes[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, absd := abs(dist_Mb)]
cand[, lab := ifelse(absd == min(absd), hit_sym, NA_integer_), by = .(qtl_id, hit_sym)]
p2c <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(data = genes, aes(dist_Mb, qtl), size = 0.35, colour = "grey65", alpha = 0.8) +
  geom_point(data = cand, aes(dist_Mb, qtl), size = 1.0, colour = "#C0392B") +
  geom_text_repel(data = cand[!is.na(lab)], aes(dist_Mb, qtl, label = lab), size = 1.9, direction = "y",
                  nudge_x = 0.3, segment.size = 0.15, min.segment.length = 0, max.overlaps = 30) +
  scale_x_continuous(limits = c(-5.2, 5.2), breaks = c(-5,-2.5,0,2.5,5)) +
  annotate("text", x = -1.4, y = "distal-ParW729-4A", label = "RBR2 / CPL3", size = 1.9, colour = "#C0392B") +
  labs(x = "Distance from peak (Mb)", y = NULL) +
  theme(axis.text.y = element_text(size = 6))
fig2 <- p2a_img / (p2b | p2c) + plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 1), widths = c(1, 1.35))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 190, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 190, units = "mm", dpi = 600)
cat("Fig2 a+b+c done\n")
