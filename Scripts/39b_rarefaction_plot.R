## 39b_rarefaction_plot.R -- raw vs rarefied CO by PD cluster
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
CO[, co := R1 + R2a + C + R2b + R3]
CO <- merge(CO, fread(file.path(TAB, "Table25_PD_clusters.csv")), by = "family")
RA <- fread(file.path(TAB, "Table28_CO_rarefied.csv"))
A <- CO[, .(cluster, co, panel = "all markers (725-3,595)")]
B <- RA[, .(cluster, co = co_rare, panel = "rarefied to 700 markers")]
D <- rbind(A, B)
D[, panel := factor(panel, levels = c("all markers (725-3,595)", "rarefied to 700 markers"))]
p <- ggplot(D, aes(x = cluster, y = co, fill = cluster)) +
  geom_violin(colour = NA, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.25) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1.2, colour = "#C0392B") +
  facet_wrap(~ panel) +
  geom_text(data = data.table(panel = factor(c("all markers (725-3,595)", "rarefied to 700 markers"),
                                             levels = levels(D$panel)),
                              label = c("P = 1.4e-09", "P = 0.83")),
            aes(x = 1.5, y = Inf, label = label), vjust = 1.6, size = 2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("#2C6E9B", "#E08A2E"), guide = "none") +
  theme_classic(base_size = 7) + labs(x = "PD cluster", y = "CO per individual")
dir.create(file.path(FIG, "supplementary"), showWarnings = FALSE)
ggsave(file.path(FIG, "supplementary", "FigS17a_CO_rarefaction.png"), p, width = 5.0, height = 2.8, dpi = 600)
cat("plot done\n")