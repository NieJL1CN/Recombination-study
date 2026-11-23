## 37_window_covariate.R -- number of informative windows per individual
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
WIN <- 15e6; STEP <- 7.5e6
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
out <- list()
for (f in COH) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  mk <- mk[order(match(chr, CHR), ppos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  nw <- integer(ncol(G))
  for (ch in unique(mk$chr)) {
    s <- mk[chr == ch]; if (nrow(s) < 5) next
    n <- ceiling((max(s$ppos) - WIN) / STEP) + 1
    if (n < 2) next
    st <- matrix(NA_real_, n, ncol(G))
    for (i in seq_len(n)) {
      a <- which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)
      if (!length(a)) next
      m <- colMeans(G[a, , drop = FALSE], na.rm = TRUE)
      st[i, ] <- ifelse(m <= 0.05, 0, ifelse(m >= 0.95, 1, NA))
    }
    nw <- nw + colSums(!is.na(st))
  }
  out[[f]] <- data.table(dataset = f, id = colnames(G), n_win = nw)
  cat("done", f, "\n")
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table23_informative_windows.csv"))
cat("rows", nrow(R), " mean n_win", round(mean(R$n_win)), "\n")
