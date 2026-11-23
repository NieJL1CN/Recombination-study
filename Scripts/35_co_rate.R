## 35_co_rate.R -- reference-style rate integration (1/span per Mb, summed per region)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- as.numeric(Sys.getenv("CO_WIN", "5")) * 1e6
STEP <- WIN / 2
SFX <- Sys.getenv("SUFFIX", "")
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
E <- merge(fread(file.path(TAB, paste0("Table18c_CO_events", SFX, ".csv"))), B, by = "chr")
## reference formula: span = Start window start -> End window end
##   = (w2 - w1) * STEP + WIN
E[, span := gap_win * STEP + WIN]
E[, start := pmax(ppos - span / 2, 0)]
E[, end := ppos + span / 2]
regs <- c("R1","R2a","C","R2b","R3")
for (r in regs) {
  lo <- switch(r, R1 = 0, R2a = E$b1, C = E$b2, R2b = E$b3, R3 = E$b4)
  hi <- switch(r, R1 = E$b1, R2a = E$b2, C = E$b3, R2b = E$b4, R3 = 1e12)
  E[, (r) := pmax(0, pmin(end, hi) - pmax(start, lo)) / span]
}
n_ind <- E[, .(n = uniqueN(id)), by = dataset]
E <- merge(E, n_ind, by = "dataset")
RG <- E[, lapply(.SD, sum), by = .(dataset, id), .SDcols = regs]
RG <- RG[, lapply(.SD, function(x) x), by = .(dataset, id)]
RG <- merge(RG, n_ind, by = "dataset")
S <- RG[, lapply(.SD, function(x) sum(x) / n[1]), by = dataset, .SDcols = regs]
S[, family := sub("_WATDE.*$", "", dataset)]
S2 <- S[, lapply(.SD, mean), by = family, .SDcols = regs]
fwrite(S2, file.path(TAB, paste0("Table22_CO_region_rate", SFX, ".csv")))
fwrite(RG, file.path(TAB, paste0("Table22c_CO_rate_individual", SFX, ".csv")))
cat("region CO per individual (rate integration):\n")
print(data.table(region = regs, co = round(sapply(S2[, ..regs], mean), 2)))
## 10 Mb bins
E[, bin_lo := floor(start / 10e6) * 10e6][, bin_hi := floor(end / 10e6) * 10e6]
BIN <- E[, .(bin = seq(bin_lo[1], bin_hi[1], by = 10e6)), by = .(dataset, id, chr, start, end, span, n)]
BIN[, w := pmax(0, pmin(bin + 10e6, end) - pmax(bin, start)) / span / n]
AGG <- BIN[, .(co = sum(w)), by = .(family = sub("_WATDE.*$", "", dataset), chr, bin)]
fwrite(AGG, file.path(TAB, paste0("Table22b_CO_10Mb_rate", SFX, ".csv")))
cat("bins:", nrow(AGG), "\n")