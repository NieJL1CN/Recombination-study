## 30_pd_tree.R -- PD tree (family names + bootstrap support) and PD heatmap
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ape)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
P <- fread(file.path(TAB, "Table19_PD_by_chr.csv"))
W <- dcast(P, family ~ chr, value.var = "pd", fill = 0)
M <- as.matrix(W[, -1]); rownames(M) <- W$family
Z <- scale(M)
h <- hclust(dist(Z), method = "average")
clade_sets <- function(hc) {
  n <- length(hc$order); cs <- vector("list", n - 1)
  for (i in seq_len(n - 1)) {
    m <- hc$merge[i, ]
    a <- if (m[1] < 0) -m[1] else cs[[m[1]]]
    b <- if (m[2] < 0) -m[2] else cs[[m[2]]]
    cs[[i]] <- sort(unique(c(a, b)))
  }
  cs
}
orig <- clade_sets(h); set.seed(1); B <- 500; boot <- integer(length(orig))
for (b in seq_len(B)) {
  idx <- sample(ncol(Z), replace = TRUE)
  cb <- clade_sets(hclust(dist(Z[, idx, drop = FALSE]), method = "average"))
  for (i in seq_along(orig)) if (any(vapply(cb, function(x) identical(x, orig[[i]]), logical(1)))) boot[i] <- boot[i] + 1
}
support <- round(100 * boot / B)
fwrite(data.table(node = seq_along(support), support = support),
       file.path(TAB, "Table24_PD_tree_bootstrap.csv"))
phy <- as.phylo(h)
png(file.path(FIG, "Fig3f_PD_tree.png"), width = 2100, height = 1700, res = 300)
plot(phy, direction = "rightwards", cex = 0.6, edge.width = 1,
     edge.col = "#2C6E9B", tip.col = "black", no.margin = FALSE,
     x.lim = c(0, max(node.depth.edgelength(phy)) * 1.35))
nodelabels(text = support, node = (length(phy$tip.label) + 1):(length(phy$tip.label) + length(support)),
           frame = "none", cex = 0.4, col = "grey25", adj = c(1.25, -0.25))
dev.off()
## heatmap ordered by the tree
H <- melt(W, id.vars = "family", variable.name = "chr", value.name = "pd")
H[, chr := factor(chr, levels = CHR)]
H[, family := factor(family, levels = W$family[h$order])]
p <- ggplot(H, aes(x = chr, y = family, fill = log10(pd + 1))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#F7F7E8", "#D9EBC0", "#9CCB8B", "#5EA97A", "#2C6E9B"), name = "log10(PD+1)") +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
ggsave(file.path(FIG, "Fig3f_PD_heatmap.png"), p, width = 5.4, height = 3.6, dpi = 600)
## two major clusters from the root split
cl <- cutree(h, k = 2)
clu <- data.table(family = names(cl), cluster = paste0("C", cl))
fwrite(clu, file.path(TAB, "Table25_PD_clusters.csv"))
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
CO[, co := R1 + R2a + C + R2b + R3]
CO <- merge(CO, clu, by = "family")
wt <- wilcox.test(co ~ cluster, data = CO)
tt <- t.test(co ~ cluster, data = CO)
cat("cluster sizes:", paste(clu[, .N, by = cluster][[2]], collapse = "/"), "\n")
cat("Wilcoxon p =", signif(wt$p.value, 3), " t-test p =", signif(tt$p.value, 3), "\n")
## ---- outlier sensitivity: drop individuals outside the Tukey fences
q <- quantile(CO$co, c(0.25, 0.75)); iqr <- q[2] - q[1]
fence <- c(q[1] - 1.5 * iqr, q[2] + 1.5 * iqr)
CO[, out := co < fence[1] | co > fence[2]]
COT <- CO[out == FALSE]
wt2 <- wilcox.test(co ~ cluster, data = COT)
tt2 <- t.test(co ~ cluster, data = COT)
cat("outliers removed:", sum(CO$out), "of", nrow(CO), "\n")
cat("trimmed Wilcoxon p =", signif(wt2$p.value, 3), " t-test p =", signif(tt2$p.value, 3), "\n")
cat("trimmed means:", paste(round(tapply(COT$co, COT$cluster, mean), 1), collapse = " / "), "\n")
fwrite(data.table(metric = c("n_all", "n_outliers", "wilcox_all", "wilcox_trimmed", "t_all", "t_trimmed",
                             "mean_C1_all", "mean_C2_all", "mean_C1_trimmed", "mean_C2_trimmed"),
                  value = c(nrow(CO), sum(CO$out), signif(wt$p.value, 3), signif(wt2$p.value, 3),
                            signif(tt$p.value, 3), signif(tt2$p.value, 3),
                            round(tapply(CO$co, CO$cluster, mean), 2),
                            round(tapply(COT$co, COT$cluster, mean), 2))),
       file.path(TAB, "Table26_cluster_CO_outlier_test.csv"))
CO[, group := "all individuals"]
COT[, group := "outliers removed"]
pl <- rbind(CO[, .(cluster, co, group)], COT[, .(cluster, co, group)])
pl[, group := factor(group, levels = c("all individuals", "outliers removed"))]
lab <- sprintf("all: P = %s\ntrimmed: P = %s", format.pval(wt$p.value, digits = 2),
               format.pval(wt2$p.value, digits = 2))
p3 <- ggplot(pl, aes(x = cluster, y = co, fill = cluster)) +
  geom_violin(colour = NA, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.25) +
  geom_jitter(data = pl, width = 0.06, size = 0.15, alpha = 0.25, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1.2, colour = "#C0392B") +
  annotate("text", x = 1.5, y = max(CO$co) * 0.98, label = lab, size = 1.9, colour = "grey15") +
  scale_fill_manual(values = c("#2C6E9B", "#E08A2E"), guide = "none") +
  facet_wrap(~ group) +
  theme_classic(base_size = 7) + labs(x = "PD cluster", y = "CO per individual")
ggsave(file.path(FIG, "Fig3f_CO_cluster_violin.png"), p3, width = 5.2, height = 3.0, dpi = 600)
print(clu[, .(families = paste(family, collapse = ", ")), by = cluster])
cat("done; support:", paste(support, collapse = ","), "\n")