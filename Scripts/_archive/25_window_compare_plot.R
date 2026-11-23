suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
D <- rbind(
  fread(file.path(TAB, "Table17_CO_direct_vs_window_20Mb.csv"))[dataset %in% COH][, win := "20 Mb / 10 Mb"],
  fread(file.path(TAB, "Table17_CO_direct_vs_window_5Mb.csv"))[dataset %in% COH][, win := "5 Mb / 2.5 Mb"])
D[, family := sub("_WATDE.*$", "", dataset)]
p <- ggplot(D, aes(x = win, y = co_window, fill = win)) +
  geom_violin(colour = NA, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.25) +
  scale_fill_manual(values = c("20 Mb / 10 Mb" = "#2C6E9B", "5 Mb / 2.5 Mb" = "#E08A2E"), guide = "none") +
  theme_classic(base_size = 8) + labs(x = NULL, y = "CO per individual")
ggsave(file.path(FIG, "Fig3h_window_5_vs_20.png"), p, width = 3.6, height = 3.4, dpi = 600)
cat("done\n")
