suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6),
  legend.title = element_blank(), plot.margin = margin(2,3,2,2)))
d <- fread(file.path(TAB, "Table172_candidate_tissue_stage.csv"))
## tissue specificity: log2(max tissue / leaf)
t <- d[set == "tissue"]
veg <- t[label %in% c("leaf","root","stem"), .(veg = mean(value)), by = gene]
t2 <- merge(t[label == "spike"], veg, by = "gene")
t2[, fold := log2((value + 0.01)/(veg + 0.01))]
tm <- t2[, .(fold = fold, tissue = "spike"), by = gene]
tm <- tm[order(-fold)]
tm[, gene := factor(gene, levels = rev(tm$gene))]
pc <- ggplot(tm, aes(fold, gene, fill = tissue)) + geom_col(width = 0.62) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c(spike = "#C0392B")) +
  labs(x = expression(log[2]~"spike / vegetative"), y = NULL) + theme(legend.position = "none")
## stage specificity: log2(max / median across stages)
s <- d[set == "stage"]
sm <- s[, .(fold = log2((max(value) + 0.01)/(median(value) + 0.01)), stage = label[which.max(value)]), by = gene]
sm <- sm[order(-fold)]
sm[, gene := factor(gene, levels = rev(sm$gene))]
pd <- ggplot(sm, aes(fold, gene, fill = stage)) + geom_col(width = 0.62) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c(latent_lepto = "#2E6FA3", zygo_pachy = "#2F8F5B", diplo_dia = "#E08214", metaphaseI = "#6A5199")) +
  labs(x = expression(log[2]~"peak stage / median"), y = NULL)
fig <- (pc | pd)
ggsave(file.path(OUT, "panels", "Fig2cd_specificity.png"), fig, width = 150, height = 62, units = "mm", dpi = 600)
cat("cd done\n")
