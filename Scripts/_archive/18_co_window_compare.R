## 18_co_window_compare.R -- CO with 5.5/2.75 vs 20/10 Mb windows
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

run_co <- function(f, WIN, STEP) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  tot <- integer(ncol(G))
  for (ch in unique(mk$chr)) {
    s <- mk[chr == ch]; if (nrow(s) < 10) next
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    W <- matrix(NA_integer_, length(starts), ncol(G))
    for (k in seq_along(starts)) {
      rr <- which(s$ppos >= starts[k] & s$ppos < starts[k] + WIN * 1e6)
      if (!length(rr)) next
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    for (j in seq_len(ncol(G))) {
      v <- W[, j]; keep <- !is.na(v); if (sum(keep) < 2) next
      kk <- which(keep); tot[j] <- tot[j] + sum(v[kk[-1]] != v[kk[-length(kk)]])
    }
  }
  data.table(dataset = f, co = as.numeric(tot))
}
res <- list()
for (w in list(c(5.5, 2.75), c(20, 10))) {
  cat("window", w[1], "step", w[2], "\n")
  for (f in COH) {
    r <- run_co(f, w[1], w[2]); if (is.null(r)) next
    r[, win := w[1]]; res[[paste(f, w[1])]] <- r
  }
}
R <- rbindlist(res)
R[, family := sub("_WATDE.*$", "", dataset)]
S <- R[, .(mean_co = round(mean(co), 1), median_co = median(co), n = .N), by = .(win, family)][order(win, -mean_co)]
fwrite(S, file.path(TAB, "Table13_CO_window_compare.csv"))
print(S[win == 5.5][1:5]); print(S[win == 20][1:5])
cat("\noverall:\n"); print(S[, .(mean = round(mean(mean_co),1), median = median(mean_co)), by = win])