suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
E <- fread(file.path(TAB, "Table178_marker_intervals.csv"))
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
n <- CO[, .(n = uniqueN(id)), by = dataset]
E <- merge(E, n, by = "dataset")
CHR <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
E <- merge(E, B, by = "chr")
regs <- c("R1","R2a","C","R2b","R3")
for (r in regs) {
  lo <- switch(r, R1 = 0, R2a = E$b1, C = E$b2, R2b = E$b3, R3 = E$b4)
  hi <- switch(r, R1 = E$b1, R2a = E$b2, C = E$b3, R2b = E$b4, R3 = 1e12)
  E[, (r) := pmax(0, pmin(b, hi) - pmax(a, lo)) / d]
}
S <- E[, lapply(.SD, sum), by = .(dataset, chr), .SDcols = regs]
S <- merge(S, n, by = "dataset")
S[, (regs) := lapply(.SD, function(x) x / n), .SDcols = regs]
F <- S[, lapply(.SD, function(x) sum(x)), by = dataset, .SDcols = regs]
F[, family := sub("_WATDE.*$", "", dataset)]
FM <- F[, lapply(.SD, mean), by = family, .SDcols = regs]
cat("region means (per RIL):\n"); print(round(sapply(FM[, ..regs], mean), 2))
cat("total per RIL:", round(sum(sapply(FM[, ..regs], mean)), 2), "\n")
# chromosome totals and correlation with length
LEN <- data.table(chr = CHR, Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
CH <- S[, lapply(.SD, function(x) sum(x)), by = .(dataset, chr), .SDcols = regs]
CH[, total := rowSums(.SD), .SDcols = regs]
CH <- merge(CH, n, by = "dataset")[, total := total / n]
CH <- CH[, .(co = mean(total)), by = chr]
CH <- merge(CH, LEN, by = "chr")
cat("r(chr length, CO):", round(cor(CH$Mb, CH$co), 3), "\n")
cat("family totals range:", round(range(FM[, rowSums(.SD), .SDcols = regs]), 2), "\n")
fwrite(FM, file.path(TAB, "Table180_region_rate_marker.csv"))
fwrite(CH, file.path(TAB, "Table181_chr_rate_marker.csv"))
