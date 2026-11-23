## 81_transmission.R -- allele transmission at significant CO-QTL peaks
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
out <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    G <- read_geno(cand); if (is.null(G)) next
    dd <- MK[dataset == cand & chr == ch]; dd <- dd[marker %in% rownames(G)]
    if (!nrow(dd)) next
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    gt <- G[pk, ]; nA <- sum(gt == 0, na.rm = TRUE); nB <- sum(gt == 1, na.rm = TRUE)
    if (nA + nB < 30) next
    high <- DIR$direction[i]
    nHigh <- if (high == "Watkins") nB else nA
    pHigh <- binom.test(nHigh, nA + nB, 0.5)$p.value
    out[[length(out) + 1]] <- data.table(family = fam, chr = ch, trait = DIR$trait[i],
      peak_marker = pk, nA_parentA = nA, nB_Watkins = nB,
      freq_A = round(nA/(nA+nB), 3), high_parent = high,
      n_high = nHigh, freq_high = round(nHigh/(nA+nB), 3), p_binom = signif(pHigh, 3))
    break
  }
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table69_QTL_allele_transmission.csv"))
print(R)