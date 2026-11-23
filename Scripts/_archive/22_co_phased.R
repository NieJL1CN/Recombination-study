## 22_co_phased.R -- phase markers against neighbours, then count CO (20/10 Mb)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
WIN <- 20e6; STEP <- 10e6

phase_chr <- function(X) {           # X: markers x individuals, 0/1
  n <- nrow(X); if (n < 3) return(X)
  ref <- X[1, ]
  for (i in 2:n) {
    a <- X[i, ]; ok <- !is.na(a) & !is.na(ref)
    if (sum(ok) < 10) { ref <- ifelse(is.na(a), ref, a); next }
    if (mean(a[ok] == ref[ok]) < 0.5) X[i, ] <- 1 - X[i, ]
    ref <- ifelse(is.na(X[i, ]), ref, X[i, ])
  }
  X
}

res <- list()
for (f in COH) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  Ph <- matrix(NA_real_, nrow(G), ncol(G))
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch); if (length(idx) < 3) { Ph[idx, ] <- G[idx, ]; next }
    Ph[idx, ] <- phase_chr(G[idx, , drop = FALSE])
  }
  tot <- integer(ncol(G))
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch); if (length(idx) < 5) next
    s <- mk[idx]; mx <- max(s$ppos); nw <- ceiling((mx - WIN) / STEP) + 1
    if (nw < 2) next
    Mx <- matrix(NA_real_, nw, ncol(G))
    for (i in seq_len(nw)) {
      a <- idx[which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)]
      if (!length(a)) next
      Mx[i, ] <- colMeans(Ph[a, , drop = FALSE], na.rm = TRUE)
    }
    Mx <- round(Mx, 1); Mx[!is.na(Mx) & Mx > 0.01 & Mx < 1.99] <- NA
    for (j in seq_len(ncol(G))) {
      v <- Mx[, j]; v <- v[!is.na(v)]; if (length(v) < 2) next
      tot[j] <- tot[j] + sum(v[-1] != v[-length(v)])
    }
  }
  res[[f]] <- data.table(dataset = f, co = as.numeric(tot))
  cat("done", f, "mean", round(mean(tot), 1), "\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table16_CO_phased.csv"))
R[, family := sub("_WATDE.*$", "", dataset)]
S <- R[, .(mean_co = mean(co), median = median(co), n = .N), by = .(dataset, family)][order(-mean_co)]
fwrite(S, file.path(TAB, "Table16b_CO_phased_summary.csv"))
print(S[1:6]); cat("overall mean:", round(mean(S$mean_co), 1), "\n")
