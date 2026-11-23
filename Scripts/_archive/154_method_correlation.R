suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds")); M <- L$markers
W <- fread(file.path(ROOT, "Tables/Table22c_CO_rate_individual.csv"))
FAMS <- c("ParW729","ParW007","ParW145","ParW352","ParW579")
out <- list()
for (fam in FAMS) {
  cand <- unique(M[grepl(paste0("^", fam, "_"), dataset), dataset])
  ds <- cand[which.max(sapply(cand, function(z) nrow(M[dataset == z & !is.na(ppos)])))]
  d <- M[dataset == ds & !is.na(ppos)]
  x <- read.csv(file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  ids <- x[[1]][-(1:2)]
  g <- x[-(1:2), -1, drop = FALSE]
  G <- t(apply(as.matrix(g), 2, function(z) { z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]
  d <- d[marker %in% rownames(G)][order(chr, ppos)]
  G <- G[d$marker, , drop = FALSE]
  co <- numeric(length(ids))
  for (ch in unique(d$chr)) {
    idx <- which(d$chr == ch); if (length(idx) < 2) next
    sub <- G[idx, , drop = FALSE]
    for (j in seq_along(ids)) { v <- sub[, j]; v <- v[!is.na(v)]; if (length(v) > 1) co[j] <- co[j] + sum(diff(v) != 0) }
  }
  w <- W[dataset == ds]
  w[, id := sub(paste0(fam, "_"), "", id)]
  id2 <- sub(paste0(fam, "_"), "", ids)
  wt <- w[match(id2, w$id)]
  win <- rowSums(wt[, .(R1, R2a, C, R2b, R3)], na.rm = TRUE)
  ok <- !is.na(win) & !is.na(co)
  out[[fam]] <- data.table(family = fam, n = sum(ok),
    pearson = round(cor(co[ok], win[ok]), 3),
    spearman = round(cor(co[ok], win[ok], method = "spearman"), 3),
    mean_trans = round(mean(co[ok]), 1), mean_window = round(mean(win[ok]), 1))
}
print(rbindlist(out))
