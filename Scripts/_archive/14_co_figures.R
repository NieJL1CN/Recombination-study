## ============================================================================
## 14_co_figures.R -- CO method comparison, violin, heatmap, genome distribution
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
REG <- c("R1" = "#2C6E9B", "R2a" = "#6A4C93", "C" = "#B0B0B0", "R2b" = "#2E8B57", "R3" = "#E08A2E")

COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
CO <- fread(file.path(TAB, "Table9_CO_by_individual.csv"))[dataset %in% COH]
CO[, family := sub("_WATDE.*$", "", dataset)]

## ---- a. method comparison bar ---------------------------------------------
M <- melt(CO[, .(marker = mean(co_marker), window = mean(co_win))], variable.name = "method", value.name = "co")
M[, method := factor(method, levels = c("marker", "window"), labels = c("marker-based", "window-based (5.5 Mb)"))]
p <- ggplot(M, aes(x = method, y = co, fill = method)) +
  geom_col(width = 0.55, colour = "grey25", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c("#2C6E9B", "#2E8B57")) +
  theme_classic(base_size = 9) + labs(x = NULL, y = "CO per individual")
ggsave(file.path(FIG, "Fig3a_CO_method_bar.png"), p, width = 4, height = 3.4, dpi = 500)

## ---- b. violin -------------------------------------------------------------
ord <- CO[, .(m = median(co_win)), by = family][order(m)]$family
V <- copy(CO); V[, family := factor(family, levels = ord)]
p2 <- ggplot(V, aes(x = family, y = co_win)) +
  geom_violin(fill = "#2E8B57", colour = NA, alpha = 0.5, scale = "width") +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.35, colour = "grey25") +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 3.2)) +
  labs(x = NULL, y = "CO per individual (window-based)")
ggsave(file.path(FIG, "Fig3b_CO_violin.png"), p2, width = 7, height = 3.6, dpi = 600)

## ---- c. chromosome x family heatmap ---------------------------------------
H <- fread(file.path(TAB, "Table10_CO_by_chr_family.csv"))
H[, chr := factor(chr, levels = CHR)]
p3 <- ggplot(H, aes(x = family, y = chr, fill = co)) + geom_tile() +
  geom_text(aes(label = sprintf("%.1f", co)), size = 1.1, colour = "grey20") +
  scale_fill_gradientn(colours = c("#F7E6A3", "#9CCB8B", "#2E8B57", "#1F3B73"), name = "mean CO") +
  theme_minimal(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 3.4),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
ggsave(file.path(FIG, "Fig3c_CO_heatmap.png"), p3, width = 9, height = 5, dpi = 600)
cat("figures a-c done\n")