## 50_co_qtl_trait_scan.R -- enrichment of CO-QTL against all published traits
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
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
WIN <- 10e6
near_frac <- function(P) {
  ## fraction of BG markers within WIN of a published QTL of this trait
  tot <- 0; hit <- 0
  for (ch in unique(BG$chr)) {
    b <- BG[chr == ch]$ppos
    p <- sort(P[chr == ch]$bp)
    if (!length(p) || !length(b)) next
    idx <- findInterval(b, p)
    left <- ifelse(idx > 0, abs(b - p[pmax(idx, 1)]), Inf)
    right <- ifelse(idx < length(p), abs(b - p[pmin(idx + 1, length(p))]), Inf)
    dmin <- pmin(left, right)
    tot <- tot + length(b); hit <- hit + sum(dmin <= WIN)
  }
  if (tot == 0) return(c(frac = NA, n = 0, k = 0))
  c(frac = hit / tot, n = tot, k = hit)
}
co_near <- function(P) {
  k <- 0
  for (i in seq_len(nrow(COQ))) {
    p <- P[chr == COQ$chr[i]]$bp
    if (length(p) && any(abs(p - COQ$phys[i]) <= WIN)) k <- k + 1
  }
  k
}
traits <- PUB[, .N, by = trait][N >= 20]$trait
res <- list()
for (tr in traits) {
  P <- PUB[trait == tr & !is.na(bp_peak), .(chr = chr2, bp = bp_peak)]
  if (!nrow(P)) next
  bg <- near_frac(P)
  k <- co_near(P); n <- nrow(COQ)
  res[[length(res) + 1]] <- data.table(trait = tr, n_pub_qtl = nrow(P),
    n_co = n, observed = k, obs_frac = k / n, expected = bg[["frac"]],
    fold = round((k / n) / max(bg[["frac"]], 1e-9), 2),
    p = pbinom(k - 1, n, bg[["frac"]], lower.tail = FALSE))
}
R <- rbindlist(res)
R[, p_adj := p.adjust(p, "BH")]
setorder(R, p)
fwrite(R, file.path(TAB, "Table39_CO_QTL_trait_scan.csv"))
cat("traits tested:", nrow(R), "\n")
print(R[1:15, .(trait, n_pub_qtl, observed, obs_frac = round(obs_frac, 3),
                expected = round(expected, 3), fold, p = signif(p, 2), p_adj = signif(p_adj, 2))])
## forest plot of fold enrichment (traits with >=30 published QTL)
lbl <- c(
  Ant_dto_day = "Anthesis / flowering date (Ant)",
  Hd_dto_day = "Heading date (Hd)",
  Mat_dto_day = "Maturity date (Mat)",
  PH_M_cm = "Plant height (PH)",
  GY_Calc_tha = "Grain yield (GY)",
  HI_Calc_pct = "Harvest index (HI)",
  BM_Calc_gm2 = "Biomass (BM)",
  NDVI_M_Idx = "NDVI",
  StrawYld_Calc_tha = "Straw yield",
  AGBMNitr_M_pct = "Above-ground N concentration",
  STRNUP_Calc_kgha = "Stem N uptake",
  BMNUtE_Calc_kgkgN = "Biomass N utilisation efficiency",
  BMNUE_Calc_kgkgN = "Biomass N use efficiency",
  NHI_Calc_pct = "N harvest index",
  NUE_Calc_kgkgN = "N use efficiency",
  NupE_Calc_kgkgN = "N uptake efficiency",
  NutE_Calc_kgkgN = "N utilisation efficiency",
  GNitr_NIRS_pct = "Grain N concentration",
  STRNitr_NIRS_pct = "Straw N concentration",
  GNUP_Calc_kgha = "Grain N uptake",
  BMNUP_Calc_kgha = "Biomass N uptake",
  GFP_calc_day = "Grain filling period",
  GFR_Calc_gday = "Grain filling rate",
  GPD_Calc_pct = "Grain protein deviation"
)
D <- R[n_pub_qtl >= 30][order(fold)]
D[, trait_label := ifelse(trait %in% names(lbl), lbl[trait], trait)]
D[, trait := factor(trait_label, levels = trait_label)]
p <- ggplot(D, aes(x = fold, y = trait)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = p_adj < 0.05, size = n_pub_qtl)) +
  scale_colour_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2C6E9B"), name = "BH < 0.05") +
  scale_size_continuous(range = c(0.4, 2), name = "published QTL") +
  scale_x_log10() + theme_classic(base_size = 7) +
  labs(title = "Colocalization with published trait QTL",
       x = "fold enrichment of CO-QTL near published QTL (log scale)", y = NULL) +
  theme(plot.title = element_text(size = 8, hjust = 0))
ggsave(file.path(FIG, "Fig4e_CO_QTL_trait_scan.png"), p, width = 4.4, height = 4.6, dpi = 600)
cat("figure done\n")