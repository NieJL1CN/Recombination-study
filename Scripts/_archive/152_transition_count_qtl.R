## 152_transition_count_qtl.R -- marker-transition counting vs window counting
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
set.seed(20261127)
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))
M <- L$markers
FAMS <- c("ParW729", "ParW007", "ParW145", "ParW352", "ParW579")
res <- list(); peaks <- list()
lod_from_r <- function(r, n) { t <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-12)); -log10(2 * pt(-abs(t), n - 2)) / 2 }
for (fam in FAMS) {
  cand_ds <- unique(M[grepl(paste0("^", fam, "_"), dataset), dataset])
  ds <- cand_ds[which.max(sapply(cand_ds, function(z) nrow(M[dataset == z & !is.na(ppos)])))]
  d <- M[dataset == ds & !is.na(ppos)]
  x <- read.csv(file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  ids <- x[[1]][-(1:2)]
  g <- x[-(1:2), -1, drop = FALSE]
  G <- t(apply(as.matrix(g), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]
  d <- d[marker %in% rownames(G)][order(chr, ppos)]
  G <- G[d$marker, , drop = FALSE]
  ## transitions per line
  co <- numeric(length(ids))
  for (ch in unique(d$chr)) {
    idx <- which(d$chr == ch); if (length(idx) < 2) next
    sub <- G[idx, , drop = FALSE]
    for (j in seq_along(ids)) {
      v <- sub[, j]; v <- v[!is.na(v)]
      if (length(v) > 1) co[j] <- co[j] + sum(diff(v) != 0)
    }
  }
  y <- co
  ## impute marker means and centre
  Gm <- G
  for (i in seq_len(nrow(Gm))) { v <- Gm[i, ]; v[is.na(v)] <- mean(v, na.rm = TRUE); Gm[i, ] <- v }
  Gm <- Gm - rowMeans(Gm)
  yc <- y - mean(y)
  denom <- sqrt(rowSums(Gm^2) * sum(yc^2))
  r <- as.vector(Gm %*% yc) / denom
  lod <- lod_from_r(r, length(y))
  maxperm <- numeric(500)
  for (b in 1:500) {
    yp <- sample(yc)
    rp <- as.vector(Gm %*% yp) / denom
    maxperm[b] <- max(lod_from_r(rp, length(y)), na.rm = TRUE)
  }
  thr <- quantile(maxperm, 0.95)
  cat(fam, "lines:", length(ids), "markers:", nrow(G), "mean CO:", round(mean(y), 1),
      "threshold:", round(thr, 2), "max LOD:", round(max(lod, na.rm=TRUE), 2), "\n")
  res[[fam]] <- data.table(family = fam, n_lines = length(ids), n_markers = nrow(G),
                           mean_co = round(mean(y), 1), threshold = round(thr, 2), max_lod = round(max(lod), 2))
  dd <- data.table(chr = d$chr, cM = d$g_rank, Mb = d$ppos/1e6, marker = d$marker, lod = lod)
  dd[, is_peak := lod == max(lod), by = chr]
  pk <- dd[is_peak == TRUE & lod >= thr]
  if (nrow(pk)) peaks[[fam]] <- data.frame(family = fam, chr = as.character(pk$chr), cM = round(pk$cM,1), Mb = round(pk$Mb,1), marker = as.character(pk$marker), lod = round(pk$lod,2), stringsAsFactors = FALSE)
}
fwrite(rbindlist(res), file.path(ROOT, "Tables/Table132_transition_qtl_summary.csv"))
if (length(peaks)) write.csv(do.call(rbind, peaks), file.path(ROOT, "Tables/Table133_transition_qtl_peaks.csv"), row.names = FALSE)
cat("done\n")
