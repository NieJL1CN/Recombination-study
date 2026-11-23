suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main"); theme_set(theme_classic(base_size = 8))
ind <- fread(file.path(TAB, "Table96_all_parw_induction.csv"))[!is.na(induction)]
ind[, hit := induction >= 1]
pc <- ggplot(ind, aes(induction, fill = hit)) + geom_histogram(binwidth = 0.25, colour = NA) +
  scale_fill_manual(values = setNames(c("grey70", "#C0392B"), c("FALSE", "TRUE")), guide = "none") +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  annotate("text", x = 1.05, y = Inf, label = paste0(sum(ind$hit), " of ", nrow(ind), " genes >= 2-fold"), vjust = 1.5, hjust = 0, size = 2.0) +
  labs(x = expression(log[2]~"meiosis / vegetative"), y = "Genes") + scale_x_continuous(limits = c(-4, 5))
ggsave(file.path(OUT, "panels", "Fig2c_induction.png"), pc, width = 85, height = 55, units = "mm", dpi = 600)
ex <- fread(file.path(TAB, "Table171_gene_stage_expr.csv"))
ex <- ex[stage %in% c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI")]
ex[, stage := factor(stage, levels = c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI"),
                     labels = c("Leptotene","Zygotene-Pachytene","Diplotene-Diakinesis","Metaphase I"))]
keep <- c("RBR2","CPL3","RECQL1","RHC1A")
ex <- ex[gene %in% keep]; ex[, gene := factor(gene, levels = keep)]
pd <- ggplot(ex, aes(stage, value, fill = gene)) + geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c(RBR2 = "#C0392B", CPL3 = "#E08214", RECQL1 = "#2E6FA3", RHC1A = "#2F8F5B")) +
  labs(x = NULL, y = expression(log[2]~"TPM")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5), legend.key.size = unit(2.5,"mm"))
ggsave(file.path(OUT, "panels", "Fig2d_stages.png"), pd, width = 95, height = 55, units = "mm", dpi = 600)
cat("panels c,d done\n")
