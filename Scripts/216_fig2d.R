suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6), legend.title = element_blank(),
  plot.margin = margin(2,3,2,2)))
ts <- fread(file.path(TAB, "Table172_candidate_tissue_stage.csv"))
keep <- c("RBR2","CPL3","RECQL1","RHC1A")
d <- ts[gene %in% keep]
d[, lab := fifelse(set == "tissue", label, label)]
d <- d[(set == "tissue" & label %in% c("leaf","spike")) | (set == "stage")]
labmap <- c(leaf = "Leaf", spike = "Spike", latent_lepto = "Leptotene", zygo_pachy = "Zygotene-Pachytene",
            diplo_dia = "Diplotene-Diakinesis", metaphaseI = "Metaphase I")
d[, lab := factor(labmap[lab], levels = unname(labmap))]
d[, gene := factor(gene, levels = keep)]
pd <- ggplot(d, aes(lab, value, fill = gene)) + geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c(RBR2 = "#C0392B", CPL3 = "#E08214", RECQL1 = "#2E6FA3", RHC1A = "#2F8F5B")) +
  labs(x = NULL, y = expression(log[2]~"TPM")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1, size = 5.2), legend.key.size = unit(2.4,"mm"))
ggsave(file.path(OUT, "panels", "Fig2d_stages.png"), pd, width = 105, height = 58, units = "mm", dpi = 600)
cat("d done\n")
