## 102_induction_all_genes.R -- expression induction for every gene in ParW QTL intervals
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
R <- fread(file.path(TAB, "Table96_all_parw_induction.csv"))
R[, qtl := paste0(family, "-", chr)]
cat("genes:", nrow(R), "  with expression:", sum(!is.na(induction)), "\n\n")
S <- R[, .(n_genes = .N, n_expr = sum(!is.na(induction)),
           n_ind1 = sum(induction >= 1, na.rm = TRUE),
           med = round(median(induction, na.rm = TRUE), 2),
           max_ind = round(max(induction, na.rm = TRUE), 2)), by = qtl][order(-n_ind1)]
print(S)
cat("\n== top 5 induced genes per QTL ==\n")
TT <- R[!is.na(induction)][order(qtl, -induction)][, head(.SD, 5), by = qtl]
for (i in seq_len(nrow(TT))) cat(sprintf("%-12s %-28s ind=%5.2f  %s\n", TT$qtl[i], TT$gene_id[i],
  TT$induction[i], substr(ifelse(is.na(TT$swiss[i]), "NA", TT$swiss[i]), 1, 45)))
R[, is_recomb := !is.na(recomb) & recomb]
p <- ggplot(R[!is.na(induction)], aes(qtl, induction)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.25) +
  geom_jitter(data = R[!is.na(induction) & is_recomb == FALSE], width = 0.18, size = 0.35, alpha = 0.35, colour = "grey55") +
  geom_point(data = R[!is.na(induction) & is_recomb == TRUE], colour = "#C0392B", size = 1.1) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = 2) +
  theme_classic(base_size = 7) + coord_flip() +
  labs(x = NULL, y = "meiosis - vegetative (log2)")
ggsave(file.path(FIG, "FigS_induction_all_genes.png"), p, width = 4.2, height = 3.0, dpi = 600)
fwrite(S, file.path(TAB, "Table97_induction_summary_by_QTL.csv"))
cat("\nfigure saved\n")
