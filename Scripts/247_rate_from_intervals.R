suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
E <- fread(file.path(TAB, "Table178_marker_intervals.csv"))
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
n <- CO[, .(n = uniqueN(id)), by = dataset]
E <- merge(E, n, by = "dataset")
W <- 10e6
E[, b0 := floor(a/W)*W][, b1 := floor(b/W)*W]
B <- E[, .(bin = seq(b0[1], b1[1], by = W)), by = .(dataset, chr, a, b, d, n)]
B[, w := pmax(0, pmin(bin + W, b) - pmax(bin, a)) / d / n]
P <- B[, .(co = sum(w)), by = .(dataset, chr, bin)][, .(co = mean(co)), by = .(dataset, chr, bin)]
A <- P[, .(co = mean(co) * 100 / 10), by = .(chr, bin)]
fwrite(A, file.path(TAB, "Table179_CO_10Mb_marker.csv"))
cat("max", round(max(A$co),2), "mean", round(mean(A$co),2), "\n")
