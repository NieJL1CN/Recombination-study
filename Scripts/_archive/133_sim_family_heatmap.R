suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
S <- fread(file.path(ROOT, "Tables/Table121_sim_family_saved.csv"))
S[, family := factor(family, levels = c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"))]
S[, dist_lab := factor(paste0(dist_cM, " cM"), levels = paste0(c(1,2,5,10,20), " cM"))]
p <- ggplot(S, aes(dist_lab, family, fill = saved)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", saved)), size = 2.2, colour = "grey15") +
  scale_fill_gradient(low = "grey92", high = "#C0392B", name = "Generations\nsaved") +
  labs(x = "Target-drag distance", y = NULL) +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "grey20"))
ggsave(file.path(ROOT, "Figures/FigS_sim_saving_heatmap.png"), p, width = 3.4, height = 2.2, dpi = 600)
cat("saved\n")
