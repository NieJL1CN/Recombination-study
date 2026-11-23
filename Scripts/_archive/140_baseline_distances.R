## 140_baseline_distances.R -- baseline fixation curves at four distances
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table123_sim_trajectories.csv"))
S <- R[family == "baseline" & gen <= 5, .(TTdd = mean(TTdd, na.rm = TRUE)), by = .(dist, gen)]
S[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(1,2,5,10), " cM"))]
p <- ggplot(S, aes(gen, TTdd, colour = dist_f)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = c("1 cM" = "#C0392B", "2 cM" = "#E08A2E",
                                 "5 cM" = "#2E8B57", "10 cM" = "#2C6E9B"), name = "Target-drag distance") +
  scale_x_continuous(breaks = 1:5, limits = c(1, 5)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_baseline_distances.png"), p,
       width = 3.6, height = 2.4, dpi = 600)
cat("done\n")
