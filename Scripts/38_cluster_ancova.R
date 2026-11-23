## 38_cluster_ancova.R -- is the PD-cluster CO difference explained by marker/window count?
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(emmeans)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
CO[, co := R1 + R2a + C + R2b + R3]
CO <- merge(CO, fread(file.path(TAB, "Table25_PD_clusters.csv")), by = "family")
CO <- merge(CO, fread(file.path(TAB, "Table23_informative_windows.csv"))[, .(dataset, id, n_win)],
            by = c("dataset", "id"))
NM <- fread(file.path(TAB, "Table0_cohort.csv"))[, .(dataset, n_marker)]
CO <- merge(CO, NM, by = "dataset")
CO[, log_nwin := log(n_win)]
CO[, log_nmark := log(n_marker)]

m0 <- lm(co ~ cluster, data = CO)
m1 <- lm(co ~ cluster + n_win, data = CO)
m2 <- lm(co ~ cluster + log_nwin, data = CO)
m3 <- lm(co ~ cluster + n_marker, data = CO)
m4 <- lm(co ~ cluster + log_nmark, data = CO)
getc <- function(m, nm = "clusterC2") {
  s <- summary(m)$coefficients
  if (!nm %in% rownames(s)) return(c(est = NA, p = NA))
  c(est = s[nm, 1], p = s[nm, 4])
}
mods <- c("co ~ cluster", "co ~ cluster + n_win", "co ~ cluster + log(n_win)",
          "co ~ cluster + n_marker", "co ~ cluster + log(n_marker)")
est <- rbind(getc(m0), getc(m1), getc(m2), getc(m3), getc(m4))
esc <- cbind(model = mods, as.data.frame(est))
names(esc) <- c("model", "cluster_effect", "p_value")
esc$cluster_effect <- round(as.numeric(as.character(esc$cluster_effect)), 2)
esc$p_value <- signif(as.numeric(as.character(esc$p_value)), 3)
esc$r2 <- round(sapply(list(m0, m1, m2, m3, m4), function(m) summary(m)$r.squared), 3)
fwrite(esc, file.path(TAB, "Table27_cluster_CO_ancova.csv"))
print(esc)
## adjusted means from the window-adjusted model
em <- emmeans(m2, ~ cluster)
print(em)
cat("\nraw means:\n"); print(CO[, .(n = .N, mean = round(mean(co), 1)), by = cluster])
