## 143_two_figures.R -- baseline and +QTL distance figures
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table127_sim_tight_distances.csv"))
S <- R[, .(TTdd = mean(TTdd, na.rm = TRUE)), by = .(family, pct, dist, gen)]
S[, dist_f := factor(paste0(dist, " cM"), levels = paste0(c(0.1,0.5,1,2,3), " cM"))]
pal <- c("0.1 cM" = "#7B241C", "0.5 cM" = "#C0392B", "1 cM" = "#E08A2E",
         "2 cM" = "#2E8B57", "3 cM" = "#2C6E9B")
base <- S[family == "baseline"]
p1 <- ggplot(base, aes(gen, TTdd, colour = dist_f)) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = pal, name = "Target-drag distance") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(7, "pt"), legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_baseline_distances.png"), p1, width = 3.8, height = 2.4, dpi = 600)
## +QTL figure: strongest QTL (ParW049, 17.3%)
q <- S[family == "ParW049"]
p2 <- ggplot() +
  geom_line(data = base, aes(gen, TTdd, colour = dist_f), linewidth = 0.3, linetype = 2) +
  geom_line(data = q, aes(gen, TTdd, colour = dist_f), linewidth = 0.45) +
  geom_point(data = q, aes(gen, TTdd, colour = dist_f), size = 0.45) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = pal, name = "Target-drag distance") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(7, "pt"), legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_QTL_distances.png"), p2, width = 3.8, height = 2.4, dpi = 600)
cat("done\n")
