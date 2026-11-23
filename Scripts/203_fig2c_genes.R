suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.5, colour = "black"),
  axis.title = element_text(size = 7.5), legend.text = element_text(size = 6),
  legend.title = element_blank(), plot.margin = margin(2,3,2,2)))
ft <- fread(file.path(TAB, "Table169_gene_features.csv"))
vr <- fread(file.path(TAB, "Table170_gene_variants.csv"))
ex <- fread(file.path(TAB, "Table171_gene_stage_expr.csv"))
ord <- c("RBR2","CPL3","RECQL1","RHC1A","ATG12")
ypos <- setNames(rev(seq_along(ord)), rev(ord))
ft[, y := ypos[gene]]; vr[, y := ypos[gene]]
body <- ft[feature != "promoter", .(lo = min(x1), hi = max(x2), y = y[1]), by = gene]
pm <- ggplot() +
  geom_segment(data = ft[feature == "promoter"], aes(x = x1, xend = x2, y = y, yend = y), colour = "grey55", linewidth = 0.4) +
  geom_segment(data = body, aes(x = lo, xend = hi, y = y, yend = y), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = ft[feature %in% c("five_prime_UTR","three_prime_UTR")], aes(xmin = x1, xmax = x2, ymin = y-0.07, ymax = y+0.07), fill = "grey55", colour = NA) +
  geom_rect(data = ft[feature == "exon"], aes(xmin = x1, xmax = x2, ymin = y-0.17, ymax = y+0.17), fill = "#2C3E50", colour = NA) +
  geom_point(data = vr, aes(x = x, y = y), colour = "#C0392B", size = 0.35, alpha = 0.6,
             position = position_jitter(height = 0.07, seed = 1)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1), labels = c("0","50","100")) +
  scale_y_continuous(breaks = ypos, labels = names(ypos), limits = c(0.4, 5.6)) +
  labs(x = "Position along locus (%)", y = NULL)
ex2 <- ex[stage %in% c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI")]
ex2[, stage := factor(stage, levels = c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI"),
                      labels = c("Leptotene","Zygotene-Pachytene","Diplotene-Diakinesis","Metaphase I"))]
ex2[, gene := factor(gene, levels = ord)]
pe <- ggplot(ex2, aes(stage, value, fill = gene)) + geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c(RBR2 = "#C0392B", CPL3 = "#E08214", RECQL1 = "#2E6FA3", RHC1A = "#2F8F5B", ATG12 = "#6A5199")) +
  labs(x = NULL, y = expression(log[2]~"TPM")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5))
fig <- (pm | pe) + plot_layout(widths = c(1.15, 1))
ggsave(file.path(OUT, "panels", "Fig2c_genes.png"), fig, width = 150, height = 62, units = "mm", dpi = 600)
cat("done\n")
