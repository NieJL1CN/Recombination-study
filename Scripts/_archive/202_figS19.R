suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main"); SUP <- file.path(ROOT, "Figures/supplementary")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.5, colour = "black"),
  axis.title = element_text(size = 7.5), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6), legend.title = element_text(size = 6), plot.margin = margin(2,3,2,2)))
## a: gene positions
genes <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
cand <- fread(file.path(TAB, "Table58_candidate_gene_report.csv"))
ord <- unique(genes[order(chr, peak_Mb), qtl_id])
genes[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, qtl := factor(qtl_id, levels = rev(ord))]
cand[, absd := abs(dist_Mb)]
cand[, lab := ifelse(absd == min(absd), hit_sym, NA_character_), by = .(qtl_id, hit_sym)]
pa <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(data = genes, aes(dist_Mb, qtl), size = 0.32, colour = "grey65", alpha = 0.8) +
  geom_point(data = cand, aes(dist_Mb, qtl), size = 0.9, colour = "#C0392B") +
  geom_text_repel(data = cand[!is.na(lab)], aes(dist_Mb, qtl, label = lab), size = 1.7, direction = "y",
                  nudge_x = 0.3, segment.size = 0.12, min.segment.length = 0, max.overlaps = 30) +
  annotate("text", x = -1.4, y = "distal-ParW729-4A", label = "RBR2 / CPL3", size = 1.8, colour = "#C0392B") +
  scale_x_continuous(limits = c(-5.2, 5.2), breaks = c(-5,-2.5,0,2.5,5)) +
  labs(x = "Distance from peak (Mb)", y = NULL, tag = "a") +
  theme(axis.text.y = element_text(size = 5.5))
## b: induction distribution
ind <- fread(file.path(TAB, "Table96_all_parw_induction.csv"))
ind <- ind[!is.na(induction)]
n63 <- nrow(ind[induction >= 1])
pb <- ggplot(ind, aes(induction)) +
  geom_histogram(binwidth = 0.25, fill = "grey55", colour = NA) +
  geom_vline(xintercept = c(0, 1), linetype = c(3, 2), linewidth = 0.35) +
  annotate("text", x = 1.05, y = Inf, label = paste0(n63, " genes >= 2-fold"), vjust = 1.5, hjust = 0, size = 2.0) +
  labs(x = expression(log[2]~"meiosis / vegetative"), y = "Genes") +
  scale_x_continuous(limits = c(-4, 5))
## c: induced genes per QTL
s <- fread(file.path(TAB, "Table97_induction_summary_by_QTL.csv"))
s[, qtl := factor(qtl, levels = s[order(n_ind1), qtl])]
pc <- ggplot(s, aes(n_ind1, qtl)) + geom_col(fill = "#4C72B0", width = 0.6) +
  geom_text(aes(label = n_ind1), hjust = -0.3, size = 1.9) +
  labs(x = "Genes induced >= 2-fold", y = NULL) + xlim(0, max(s$n_ind1) * 1.2) +
  theme(axis.text.y = element_text(size = 5.5))
fig <- (pa / (pb | pc)) + plot_annotation(tag_levels = list(c("a","b","c"))) + plot_layout(heights = c(1.1, 1))
ggsave(file.path(SUP, "FigS19_candidate_expression.png"), fig, width = 160, height = 150, units = "mm", dpi = 600)
cat("S19 done; n63 =", n63, "\n")
