suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
S <- fread(file.path(ROOT, "Tables/Table118_sim_replicates_summary.csv"))
base <- S[mult == 1, .(base_cM, base_gen = gen_TTdd_mean)]
S <- merge(S, base, by = "base_cM")
S[, saved := base_gen - gen_TTdd_mean]
S[, mult_lab := factor(mult, levels = c(1, 1.02, 1.17, 2, 5, 10),
                       labels = c("x1", "x1.02", "x1.17", "x2", "x5", "x10"))]
S[, dist_lab := factor(paste0(base_cM, " cM"), levels = paste0(c(1,2,5,10,20,50), " cM"))]
p <- ggplot(S, aes(mult_lab, dist_lab, fill = saved)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", saved)), size = 2.1, colour = "grey15") +
  scale_fill_gradient2(low = "grey92", mid = "grey80", high = "#C0392B", midpoint = 0.5,
                       name = "Generations\nsaved") +
  labs(x = "Increase in recombination", y = "Target-drag distance") +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "grey20"))
ggsave(file.path(ROOT, "Figures/FigS_sim_saving_heatmap.png"), p, width = 3.2, height = 2.4, dpi = 600)
cat("heatmap saved\n")
