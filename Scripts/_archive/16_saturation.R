## 16_saturation.R -- markers vs CO, saturating fit, 500-marker threshold
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CO <- fread(file.path(TAB, "Table9_CO_by_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
S <- CO[, .(co = mean(co_win), co_marker = mean(co_marker)), by = dataset]
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
n <- M[, .(n_marker = .N), by = dataset]
S <- merge(S, n, by = "dataset")
S[, cohort := ifelse(n_marker >= 500, ">=500 markers", "<500 markers")]
f1 <- try(nls(co ~ a - b * exp(-c * n_marker), data = S,
              start = list(a = 300, b = 280, c = 0.002)), silent = TRUE)
plateau <- NA_real_; r2 <- NA_real_
if (!inherits(f1, "try-error")) {
  plateau <- coef(f1)["a"]; r2 <- 1 - sum(residuals(f1)^2) / sum((S$co - mean(S$co))^2)
}
p <- ggplot(S, aes(x = n_marker, y = co, colour = cohort)) +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_smooth(method = "nls", formula = y ~ a - b * exp(-c * x),
              method.args = list(start = list(a = 300, b = 280, c = 0.002)),
              se = FALSE, colour = "grey25", linewidth = 0.6) +
  geom_vline(xintercept = 500, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  scale_colour_manual(values = c(">=500 markers" = "#2C6E9B", "<500 markers" = "#E08A2E"), name = NULL) +
  theme_classic(base_size = 8) +
  labs(x = "markers per family", y = "CO per individual",
       subtitle = sprintf("plateau = %.0f CO, R2 = %.2f", plateau, r2))
dir.create(file.path(FIG, "supplementary"), showWarnings = FALSE)
ggsave(file.path(FIG, "supplementary", "FigS17b_marker_CO_saturation.png"), p, width = 4.6, height = 3.6, dpi = 600)
fwrite(S, file.path(TAB, "Table12_markers_vs_CO.csv"))
print(S[, .(families = .N, mean_co = round(mean(co))), by = cohort])
cat("plateau", round(plateau), "r2", round(r2, 2), "\n")
