## 42_co_qtl_fitness.R -- within-family allele effects of CO-QTL on fitness traits
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, g = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
PH <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
PK[, family := sub("_WATDE.*$", "", family)]
TRAITS <- c("PH_M_cm", "GY_Calc_tha", "GW_M_g1000grn")
cache <- list(); out <- list()
for (i in seq_len(nrow(PK))) {
  fam <- PK$family[i]; ch <- PK$chr[i]; pos <- PK$pos[i]
  ds <- unique(MK[family == fam & chr == ch, dataset])
  if (!length(ds)) next
  ds <- ds[1]
  if (is.null(cache[[ds]])) {
    x <- read.csv(file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
    g <- x[-(1:2), ]
    G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
      z <- trimws(z); o <- rep(NA_real_, length(z))
      o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
    rownames(G) <- colnames(x)[-1]
    cache[[ds]] <- list(G = G, ids = g[[1]])
  }
  G <- cache[[ds]]$G; ids <- cache[[ds]]$ids
  ## nearest marker to the CO-QTL peak on this family map (same dataset)
  cand <- MK[dataset == ds & chr == ch]
  if (!nrow(cand)) next
  best <- cand[which.min(abs(g - pos))]
  if (!best$marker %in% rownames(G)) next
  gt <- G[best$marker, ]
  ph <- PH[line %in% ids]
  if (!nrow(ph)) next
  ph <- ph[match(ids, line)]
  for (tr in TRAITS) {
    y <- ph[[tr]]
    if (sum(!is.na(y) & !is.na(gt)) < 40) next
    a <- y[!is.na(gt) & gt == 0]; b <- y[!is.na(gt) & gt == 1]
    if (length(a) < 15 || length(b) < 15) next
    tt <- t.test(b, a)
    out[[length(out) + 1]] <- data.table(
      family = fam, chr = ch, co_trait = PK$trait[i], pos_cM = round(pos, 1),
      marker = best$marker, trait = tr, n0 = length(a), n1 = length(b),
      mean_paragon = round(mean(a, na.rm = TRUE), 3), mean_watkins = round(mean(b, na.rm = TRUE), 3),
      effect = round(mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE), 3),
      p = tt$p.value)
  }
}
R <- rbindlist(out)
R[, p_adj := p.adjust(p, "BH")]
fwrite(R, file.path(TAB, "Table31_CO_QTL_fitness_allele_effects.csv"))
cat("tests:", nrow(R), " significant (BH<0.05):", sum(R$p_adj < 0.05), "\n")
print(R[p_adj < 0.05][order(p_adj)][, .(family, chr, co_trait, trait, effect, p_adj = signif(p_adj, 2))])