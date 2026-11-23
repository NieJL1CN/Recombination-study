suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
S <- fread(file.path(TAB, "Table13_CO_window_compare.csv"))
S[, winlab := factor(win, levels = c(5.5, 20), labels = c("5.5 Mb window / 2.75 Mb step", "20 Mb window / 10 Mb step"))]
p <- ggplot(S, aes(x = winlab, y = mean_co)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, colour = "grey35", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.2, colour = "#2C6E9B", alpha = 0.85) +
  theme_classic(base_size = 8) +
  labs(x = NULL, y = "mean CO per individual (family)")
ggsave(file.path(FIG, "Fig3g_window_size_compare.png"), p, width = 4.4, height = 3.2, dpi = 600)
cat("done\n")
