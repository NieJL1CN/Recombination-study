## 24_co_plots.R -- method comparison, saturation, violin (corrected data)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
IND <- fread(file.path(TAB, "Table17_CO_direct_vs_window.csv"))
S <- fread(file.path(TAB, "Table17b_CO_family_summary.csv"))
COV <- fread(file.path(TAB, "Table0_cohort.csv"))[, .(dataset, n_marker)]
if (!"n_marker" %in% names(S)) S <- merge(S, COV, by = "dataset", all.x = TRUE)
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
S[, cohort := dataset %in% COH]
IND[, cohort := dataset %in% COH]

## a. method comparison bar (cohort)
M <- S[cohort == TRUE, .(marker = mean(co_marker), window = mean(co_window))]
M <- melt(M, variable.name = "method", value.name = "co")
p <- ggplot(M, aes(x = method, y = co, fill = method)) +
  geom_col(width = 0.5, colour = "grey30", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c("#2C6E9B", "#2E8B57")) +
  theme_classic(base_size = 8) + labs(x = NULL, y = "CO per individual")
ggsave(file.path(FIG, "Fig3a_CO_method_bar.png"), p, width = 3.4, height = 3.2, dpi = 600)

## b. violin (cohort, window method)
V <- IND[cohort == TRUE]; V[, family := sub("_WATDE.*$", "", dataset)]
ord <- V[, .(m = median(co_window)), by = family][order(m)]$family
V[, family := factor(family, levels = ord)]
p2 <- ggplot(V, aes(x = family, y = co_window)) +
  geom_violin(fill = "#2E8B57", colour = NA, alpha = 0.5, scale = "width") +
  geom_jitter(width = 0.12, size = 0.25, alpha = 0.5, colour = "grey25") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5)) +
  labs(x = NULL, y = "CO per individual (window method)")
ggsave(file.path(FIG, "Fig3b_CO_violin.png"), p2, width = 5.6, height = 3.6, dpi = 600)

## c. saturation: CO vs markers
f1 <- try(nls(co_window ~ a - b * exp(-c * n_marker), data = S,
              start = list(a = 50, b = 45, c = 0.002)), silent = TRUE)
plateau <- NA; r2 <- NA
if (!inherits(f1, "try-error")) {
  plateau <- coef(f1)["a"]; r2 <- 1 - sum(residuals(f1)^2) / sum((S$co_window - mean(S$co_window))^2)
}
p3 <- ggplot(S, aes(x = n_marker, y = co_window)) +
  geom_point(aes(colour = cohort), size = 1.1, alpha = 0.8) +
  geom_smooth(method = "nls", formula = y ~ a - b * exp(-c * x),
              method.args = list(start = list(a = 50, b = 45, c = 0.002)),
              se = FALSE, colour = "grey25", linewidth = 0.6) +
  geom_vline(xintercept = 500, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  scale_colour_manual(values = c("TRUE" = "#2C6E9B", "FALSE" = "#E08A2E"),
                      labels = c("TRUE" = ">=500 markers", "FALSE" = "<500 markers"), name = NULL) +
  theme_classic(base_size = 8) +
  labs(x = "markers per family", y = "CO per individual",
       subtitle = sprintf("plateau = %.0f CO, R2 = %.2f", plateau, r2))
dir.create(file.path(FIG, "supplementary"), showWarnings = FALSE)
ggsave(file.path(FIG, "supplementary", "FigS17b_marker_CO_saturation.png"), p3, width = 4.6, height = 3.6, dpi = 600)
fwrite(S, file.path(TAB, "Table17b_CO_family_summary.csv"))
cat("plateau", round(plateau), "r2", round(r2,2), "\n")
print(S[cohort == TRUE, .(n = .N, marker = round(mean(co_marker)), window = round(mean(co_window)))])