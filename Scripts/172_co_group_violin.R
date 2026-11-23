suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
D <- fread(file.path(ROOT, "Tables/Table162_co_by_line_groups.csv"))
P <- D[, .(mean = mean(co), n = .N), by = .(group, pop)]
cat("pop means sd:", P[, .(sd = round(sd(mean),2), n = .N), by = group]$sd, "\n")
cat("t test:", t.test(mean ~ group, data = P)$p.value, "\n")
cat("var test:", var.test(mean ~ group, data = P)$p.value, "\n")
cat("wilcox lines:", wilcox.test(co ~ group, data = D)$p.value, "\n")
p <- ggplot(D, aes(x = group, y = co, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Crossovers per RIL") +
  scale_fill_manual(values = c(landrace = "#4C72B0", cultivar = "#DD8452")) +
  theme(legend.position = "none")
ggsave(file.path(ROOT, "Figures/FigS18_CO_landrace_vs_cultivar.png"), p, width = 4.2, height = 4.2, dpi = 600)
cat("saved\n")
