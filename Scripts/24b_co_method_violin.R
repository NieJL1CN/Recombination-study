## 24b_co_method_violin.R -- Fig3a: distributions of CO by counting method
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
D <- fread(file.path(TAB, "Table17_CO_direct_vs_window_20Mb.csv"))[dataset %in% COH]
L <- melt(D[, .(dataset, id, marker = co_marker, window = co_window)],
          id.vars = c("dataset", "id"), variable.name = "method", value.name = "co")
L[, method := factor(method, levels = c("marker", "window"),
                     labels = c("direct marker count", "sliding window (5 Mb)"))]
p <- ggplot(L, aes(x = method, y = co, fill = method)) +
  geom_violin(colour = NA, alpha = 0.5, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.25) +
  geom_jitter(width = 0.05, size = 0.12, alpha = 0.25, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1.1, colour = "#C0392B") +
  scale_fill_manual(values = c("direct marker count" = "#2C6E9B", "sliding window (5 Mb)" = "#2E8B57"), guide = "none") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(size = 6)) +
  labs(x = NULL, y = "CO per individual")
ggsave(file.path(FIG, "Fig3a_CO_method_bar.png"), p, width = 3.2, height = 3.2, dpi = 600)
cat("Fig3a rebuilt as violin; n =", nrow(L), "\n")
