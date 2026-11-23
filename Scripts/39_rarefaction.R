## 39_rarefaction.R -- recompute CO with an equal number of markers (700)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
WIN <- as.numeric(Sys.getenv("CO_WIN", "15")) * 1e6; STEP <- WIN / 2; NSUB <- 700
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
POS <- unique(M[, .(marker, chr, ppos)])
set.seed(42)
out <- list()
for (f in COH) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), POS, by = "marker")
  mk <- mk[order(match(chr, CHR), ppos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  if (nrow(mk) > NSUB) {
    keep <- sort(sample(nrow(mk), NSUB))
    mk <- mk[keep]; G <- G[keep, , drop = FALSE]
  }
  n <- ncol(G); tot <- numeric(n)
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch); if (length(idx) < 5) next
    s <- mk[idx]; mx <- max(s$ppos); nw <- ceiling((mx - WIN) / STEP) + 1
    if (nw < 2) next
    st <- matrix(NA_real_, nw, n)
    for (i in seq_len(nw)) {
      a <- which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)
      if (!length(a)) next
      m <- colMeans(G[idx[a], , drop = FALSE], na.rm = TRUE)
      st[i, ] <- ifelse(m <= 0.05, 0, ifelse(m >= 0.95, 1, NA))
    }
    for (j in seq_len(n)) {
      v <- st[, j]; v <- v[!is.na(v)]
      if (length(v) >= 2) tot[j] <- tot[j] + sum(v[-1] != v[-length(v)])
    }
  }
  out[[f]] <- data.table(dataset = f, id = colnames(G), co_rare = tot)
  cat("done", f, "markers used", nrow(mk), "\n")
}
R <- rbindlist(out)
R[, family := sub("_WATDE.*$", "", dataset)]
R <- merge(R, fread(file.path(TAB, "Table25_PD_clusters.csv")), by = "family")
wt <- wilcox.test(co_rare ~ cluster, data = R)
cat("\nrarefied (700 markers): Wilcoxon p =", signif(wt$p.value, 3), "\n")
print(R[, .(n = .N, mean = round(mean(co_rare), 1), median = median(co_rare)), by = cluster])
fwrite(R, file.path(TAB, "Table28_CO_rarefied.csv"))