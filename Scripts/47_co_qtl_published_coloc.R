## 47_co_qtl_published_coloc.R -- CO-QTL vs published trait QTL (binomial enrichment)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
COQ <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
PUB <- fread(file.path(TAB, "Table35_published_QTL_cohort.csv"))
PUB[, chr2 := sub("_.*$", "", chr)]
COQ[, phys := NA_real_]
for (i in seq_len(nrow(COQ))) {
  d <- MK[family == COQ$family[i] & chr == COQ$chr[i]]
  if (nrow(d)) COQ$phys[i] <- d[which.min(abs(cM - COQ$pos[i])), ppos]
}
COQ <- COQ[!is.na(phys)]
BG <- unique(MK[, .(chr, ppos)])
window <- 10e6
res <- list()
for (tr in c("Ant_dto_day","Hd_dto_day","PH_M_cm","GY_Calc_tha","Mat_dto_day",
             "GNitr_NIRS_pct","NUE_Calc_kgkgN","HI_Calc_pct")) {
  P <- PUB[trait == tr & !is.na(bp_peak), .(chr = chr2, bp = bp_peak)]
  if (!nrow(P)) next
  bg_hit <- BG[, hit := sapply(seq_len(.N), function(i) any(abs(P[chr == BG$chr[i]]$bp - ppos[i]) <= window))]
  p_bg <- mean(bg_hit$hit)
  co_hit <- COQ[, hit := sapply(seq_len(.N), function(i) any(abs(P[chr == COQ$chr[i]]$bp - phys[i]) <= window))]
  k <- sum(co_hit$hit); n <- nrow(co_hit)
  res[[length(res) + 1]] <- data.table(
    trait = tr, n_co_qtl = n, observed = k, obs_frac = round(k / n, 3),
    expected_frac = round(p_bg, 3), fold = round((k / n) / max(p_bg, 1e-9), 2),
    p_binom = pbinom(k - 1, n, p_bg, lower.tail = FALSE))
}
R <- rbindlist(res)
R[, p_adj := p.adjust(p_binom, "BH")]
fwrite(R, file.path(TAB, "Table36_CO_QTL_coloc_published.csv"))
print(R[order(p_binom)])
