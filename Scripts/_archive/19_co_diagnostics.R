## 19_co_diagnostics.R -- why are CO counts high? segment length / singleton windows
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- 20; STEP <- 10
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
fams <- c("ParW471_WATDE0062_watseq_axiom", "ParW034_WATDE0005_watseq_axiom", "ParW729_WATDE0100_watseq_axiom")
for (f in fams) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  nwin <- 0; nmark <- 0; sing <- 0; seg <- list(); co <- integer(ncol(G))
  for (ch in unique(mk$chr)) {
    s <- mk[chr == ch]; if (nrow(s) < 10) next
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    W <- matrix(NA_integer_, length(starts), ncol(G)); kk <- integer(length(starts))
    for (k in seq_along(starts)) {
      rr <- which(s$ppos >= starts[k] & s$ppos < starts[k] + WIN * 1e6)
      kk[k] <- length(rr); nmark <- nmark + length(rr); nwin <- nwin + 1
      if (!length(rr)) next
      if (length(rr) == 1) sing <- sing + 1
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    for (j in seq_len(ncol(G))) {
      v <- W[, j]; keep <- !is.na(v); if (sum(keep) < 2) next
      id <- which(keep); sw <- id[-1][v[id[-1]] != v[id[-length(id)]]]
      co[j] <- co[j] + length(sw)
      if (length(sw)) seg[[length(seg) + 1]] <- data.table(family = f, chr = ch,
        seg_Mb = diff(c(starts[sw], max(s$ppos) / 1e6)))
    }
  }
  cat("\n==", f, " mean CO", round(mean(co),1), " windows", nwin,
      " singleton-window fraction", round(sing / nwin, 3),
      " markers/window", round(nmark / nwin, 1), "\n")
  S <- rbindlist(seg)
  cat("  segment length quantiles (Mb):", paste(round(quantile(S$seg_Mb, c(.1,.25,.5,.75,.9)),1), collapse=" | "), "\n")
  cat("  segments < 20 Mb:", round(mean(S$seg_Mb < 20), 3), " < 10 Mb:", round(mean(S$seg_Mb < 10), 3), "\n")
}
