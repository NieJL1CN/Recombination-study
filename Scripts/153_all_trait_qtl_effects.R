## 153_all_trait_qtl_effects.R
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, g = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
PH <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
PK <- fread(file.path(TAB, "Table20b_CO_QTL_peaks.csv"))
PK[, family := sub("_WATDE.*$", "", dataset)]
TRAITS <- c("PH_M_cm","GY_Calc_tha","GW_M_g1000grn","GN_Calc_Grnm2","GWid_M_mm","GLng_M_mm")
cache <- list(); out <- list()
for (i in seq_len(nrow(PK))) {
  fam <- PK$family[i]; ch <- PK$chr[i]; pos <- PK$pos[i]
  ds <- unique(MK[family == fam & chr == ch, dataset]); if (!length(ds)) next
  ds <- ds[1]
  if (is.null(cache[[ds]])) {
    x <- read.csv(file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
    g <- x[-(1:2), ]
    G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
      z <- trimws(z); o <- rep(NA_real_, length(z))
      o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
    rownames(G) <- colnames(x)[-1]; cache[[ds]] <- list(G = G, ids = g[[1]])
  }
  G <- cache[[ds]]$G; ids <- cache[[ds]]$ids
  cand <- MK[dataset == ds & chr == ch]; if (!nrow(cand)) next
  best <- cand[which.min(abs(g - pos))]
  if (!best$marker %in% rownames(G)) next
  gt <- G[best$marker, ]
  ph <- PH[line %in% ids]; if (!nrow(ph)) next
  ph <- ph[match(ids, line)]
  for (tr in TRAITS) {
    y <- ph[[tr]]; ok <- !is.na(y) & !is.na(gt)
    a <- y[ok & gt == 0]; b <- y[ok & gt == 1]
    if (length(a) < 15 || length(b) < 15) next
    tt <- t.test(b, a)
    out[[length(out) + 1]] <- data.table(qtl = paste0(fam,"-",ch), family = fam, chr = ch,
      co_trait = PK$trait[i], trait = tr, n_paragon = length(a), n_watkins = length(b),
      mean_paragon = round(mean(a),3), mean_watkins = round(mean(b),3),
      effect = round(mean(b) - mean(a),3), p = tt$p.value)
  }
}
R <- rbindlist(out)
R[, p_adj := p.adjust(p, "BH")]
fwrite(R, file.path(TAB, "Table134_all_trait_allele_effects.csv"))
cat("tests:", nrow(R), " significant (BH<0.05):", sum(R$p_adj < 0.05), "\n")
print(R[order(p)][1:6, .(qtl, trait, effect, p = signif(p,2), p_adj = signif(p_adj,2))])
