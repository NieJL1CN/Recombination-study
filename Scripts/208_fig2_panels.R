suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6), legend.title = element_blank(),
  plot.margin = margin(2,3,2,2), plot.tag = element_text(size = 9, face = "bold")))
d <- fread(file.path(TAB, "Table172_candidate_tissue_stage.csv"))
## tissue
t <- d[set == "tissue"]
veg <- t[label %in% c("leaf","root","stem"), .(veg = mean(value)), by = gene]
t2 <- merge(t[label == "spike"], veg, by = "gene"); t2[, fold := log2((value + 0.01)/(veg + 0.01))]
tm <- t2[order(-fold)]; tm[, gene := factor(gene, levels = rev(tm$gene))]
pc <- ggplot(tm, aes(fold, gene)) + geom_col(fill = "#C0392B", width = 0.6) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  labs(x = expression(log[2]~"spike / vegetative"), y = NULL)
## stage
s <- d[set == "stage"]
sm <- s[, .(fold = log2((max(value)+0.01)/(median(value)+0.01)), stage = label[which.max(value)]), by = gene]
sm <- sm[order(-fold)]; sm[, gene := factor(gene, levels = rev(sm$gene))]
pd <- ggplot(sm, aes(fold, gene, fill = stage)) + geom_col(width = 0.6) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c(latent_lepto = "#2E6FA3", zygo_pachy = "#2F8F5B", diplo_dia = "#E08214", metaphaseI = "#6A5199")) +
  labs(x = expression(log[2]~"peak stage / median"), y = NULL) + theme(legend.position = "right")
ggsave(file.path(OUT, "panels", "Fig2c_tissue.png"), pc, width = 78, height = 62, units = "mm", dpi = 600)
ggsave(file.path(OUT, "panels", "Fig2d_stage.png"), pd, width = 95, height = 62, units = "mm", dpi = 600)
## gene models only
ft <- fread(file.path(TAB, "Table169_gene_features.csv")); vr <- fread(file.path(TAB, "Table170_gene_variants.csv"))
ord <- c("RBR2","CPL3","RECQL1","RHC1A","ATG12"); ypos <- setNames(rev(seq_along(ord)), rev(ord))
ft[, y := ypos[gene]]; vr[, y := ypos[gene]]
body <- ft[feature != "promoter", .(lo = min(x1), hi = max(x2), y = y[1]), by = gene]
pe <- ggplot() +
  geom_segment(data = ft[feature == "promoter"], aes(x = x1, xend = x2, y = y, yend = y), colour = "grey55", linewidth = 0.4) +
  geom_segment(data = body, aes(x = lo, xend = hi, y = y, yend = y), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = ft[feature %in% c("five_prime_UTR","three_prime_UTR")], aes(xmin = x1, xmax = x2, ymin = y-0.07, ymax = y+0.07), fill = "grey55", colour = NA) +
  geom_rect(data = ft[feature == "exon"], aes(xmin = x1, xmax = x2, ymin = y-0.17, ymax = y+0.17), fill = "#2C3E50", colour = NA) +
  geom_point(data = vr, aes(x = x, y = y), colour = "#C0392B", size = 0.4, alpha = 0.6, position = position_jitter(height = 0.07, seed = 1)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1), labels = c("0","50","100")) +
  scale_y_continuous(breaks = ypos, labels = names(ypos), limits = c(0.4, 5.6)) +
  labs(x = "Position along locus (%)", y = NULL) +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(OUT, "panels", "Fig2e_models.png"), pe, width = 150, height = 55, units = "mm", dpi = 600)
cat("panels done\n")
