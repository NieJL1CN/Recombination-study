## 147_breakage_combined.R -- linkage breakage across 0.1-3 cM, baseline vs QTL
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table127_sim_tight_distances.csv"))
S <- R[family %in% c("baseline", "ParW049"),
       .(cleanT = mean(cleanT, na.rm = TRUE)), by = .(family, dist, gen)]
S[, type := fifelse(family == "baseline", "No QTL", "With QTL (ParW049, +17.3%)")]
S[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(0.1,0.5,1,2,3), " cM"))]
pal <- c("0.1 cM" = "#7B241C", "0.5 cM" = "#C0392B", "1 cM" = "#E08A2E",
         "2 cM" = "#2E8B57", "3 cM" = "#2C6E9B")
p <- ggplot(S, aes(gen, cleanT, colour = dist_f, linetype = type, group = interaction(dist_f, type))) +
  geom_line(linewidth = 0.42) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = pal, name = "Target-drag distance") +
  scale_linetype_manual(values = c("No QTL" = "solid", "With QTL (ParW049, +17.3%)" = "22"), name = NULL) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "T carriers free of the deleterious allele") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 5.4), legend.key.width = unit(14, "pt"))
ggsave(file.path(ROOT, "Figures/FigS_breakage_QTL_vs_baseline.png"), p, width = 4.0, height = 2.5, dpi = 600)
cat("done\n")
