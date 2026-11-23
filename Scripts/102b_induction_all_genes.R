## 102b_induction_all_genes.R -- induction summary using base R
suppressPackageStartupMessages({library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
R <- read.csv(file.path(TAB, "Table96_all_parw_induction.csv"), stringsAsFactors = FALSE)
R <- R[!is.na(R$induction), ]
R$qtl <- paste0(R$family, "-", R$chr)
cat("genes with expression:", nrow(R), "\n\n")
agg <- aggregate(induction ~ qtl, R, function(x) c(n = length(x), ind1 = sum(x >= 1), med = median(x), mx = max(x)))
S <- data.frame(qtl = agg$qtl, n = agg$induction[, "n"], n_ind1 = agg$induction[, "ind1"],
                med = round(agg$induction[, "med"], 2), max = round(agg$induction[, "mx"], 2))
S <- S[order(-S$n_ind1), ]
print(S, row.names = FALSE)
cat("\n== top 5 induced per QTL ==\n")
for (q in unique(R$qtl)) {
  d <- R[R$qtl == q, ]
  d <- d[order(-d$induction), ][1:min(5, nrow(d)), ]
  for (i in seq_len(nrow(d))) cat(sprintf("%-12s %-28s ind=%5.2f  %s\n", q, d$gene_id[i],
    d$induction[i], substr(ifelse(is.na(d$swiss[i]), "NA", d$swiss[i]), 1, 45)))
}
R$is_recomb <- !is.na(R$recomb) & R$recomb
p <- ggplot(R, aes(qtl, induction)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.25) +
  geom_jitter(data = R[!R$is_recomb, ], width = 0.18, size = 0.35, alpha = 0.35, colour = "grey55") +
  geom_point(data = R[R$is_recomb, ], colour = "#C0392B", size = 1.1) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = 2) +
  theme_classic(base_size = 7) + coord_flip() +
  labs(x = NULL, y = "meiosis - vegetative (log2)")
ggsave(file.path(FIG, "FigS_induction_all_genes.png"), p, width = 4.2, height = 3.0, dpi = 600)
write.csv(S, file.path(TAB, "Table97_induction_summary_by_QTL.csv"), row.names = FALSE)
cat("\nfigure saved\n")
