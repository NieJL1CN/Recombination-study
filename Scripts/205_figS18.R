suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); SUP <- file.path(ROOT, "Figures/supplementary")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.5, colour = "black"),
  axis.title = element_text(size = 7.5), plot.margin = margin(2,3,2,2)))
genes <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
cand <- fread(file.path(TAB, "Table58_candidate_gene_report.csv"))
ord <- unique(genes[order(chr, peak_Mb), qtl_id])
genes[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, absd := abs(dist_Mb)]
cand[, lab := ifelse(absd == min(absd), hit_sym, NA_character_), by = .(qtl_id, hit_sym)]
p <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(data = genes, aes(dist_Mb, qtl), size = 0.35, colour = "grey65", alpha = 0.8) +
  geom_point(data = cand, aes(dist_Mb, qtl), size = 1.0, colour = "#C0392B") +
  geom_text_repel(data = cand[!is.na(lab)], aes(dist_Mb, qtl, label = lab), size = 1.8, direction = "y",
                  nudge_x = 0.3, segment.size = 0.12, min.segment.length = 0, max.overlaps = 30) +
  annotate("point", x = -0.23, y = "distal-ParW729-4A", colour = "#C0392B", size = 1.4) +
  annotate("text", x = -0.23, y = "distal-ParW729-4A", label = "RBR2", size = 1.8, colour = "#C0392B", hjust = -0.2, vjust = -1.2) +
  annotate("point", x = 1.32, y = "distal-ParW729-4A", colour = "#C0392B", size = 1.4) +
  annotate("text", x = 1.32, y = "distal-ParW729-4A", label = "CPL3", size = 1.8, colour = "#C0392B", hjust = -0.2, vjust = -1.2) +
  scale_x_continuous(limits = c(-5.2, 5.2), breaks = c(-5,-2.5,0,2.5,5)) +
  labs(x = "Distance from peak (Mb)", y = NULL) +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(SUP, "FigS18_interval_genes.png"), p, width = 160, height = 105, units = "mm", dpi = 600)
cat("S18 done\n")
