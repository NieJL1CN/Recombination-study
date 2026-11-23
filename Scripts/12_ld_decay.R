## ============================================================================
## 12_ld_decay.R -- genome-wide LD decay, window size from r2 = 0.8
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))

L <- readRDS(file.path(DER, "maps_min1.rds")); M <- L$markers
M <- M[!is.na(g_rank) & panel == "AXIOM.1"]
pos <- unique(M[, .(marker, chr, ppos)])

fs <- unique(M$dataset)[1:6]
res <- list()
for (f in fs) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv"))
  if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]
  mk <- data.table(marker = rownames(G))
  mk <- merge(mk, pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  for (ch in unique(mk$chr)) {
    sel <- which(mk$chr == ch); if (length(sel) < 30) next
    sel <- sel[order(mk$ppos[sel])]
    sel <- sel[seq(1, length(sel), by = max(1, floor(length(sel) / 400)))]
    X <- G[sel, , drop = FALSE]
    if (nrow(X) < 20) next
    P <- mk$ppos[sel]
    C <- cor(t(X), use = "pairwise.complete.obs")
    d <- abs(outer(P, P, "-")) / 1e6
    iu <- upper.tri(C)
    res[[length(res) + 1]] <- data.table(family = f, chr = ch,
      dist_Mb = d[iu], r2 = as.numeric(C[iu])^2)
  }
  cat("done", f, "\n")
}
D <- rbindlist(res)
D <- D[dist_Mb > 0 & dist_Mb < 100]
D[, bin := cut(dist_Mb, breaks = seq(0, 100, by = 1), labels = FALSE)]
S <- D[, .(r2 = mean(r2, na.rm = TRUE), n = .N), by = bin][order(bin)]
S[, dist := bin - 0.5]
fwrite(S, file.path(TAB, "Table7_LD_decay.csv"))
cat("\nLD by distance (first bins):\n"); print(head(S, 10))
w80 <- S[r2 <= 0.8][1, dist]
cat("window where r2 <= 0.8:", w80, "Mb\n")
cat("half-decay distance:", S[r2 <= S[1, r2] / 2][1, dist], "Mb\n")
fwrite(data.table(window_Mb = w80, half_Mb = S[r2 <= S[1, r2] / 2][1, dist]),
       file.path(TAB, "Table8_window_from_LD.csv"))