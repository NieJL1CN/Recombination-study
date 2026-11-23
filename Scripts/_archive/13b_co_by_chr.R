## 13b_co_by_chr_family.R -- mean window CO per family x chromosome
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- 5.5; STEP <- 2.75
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers
M <- M[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
out <- list()
for (f in unique(M$dataset)) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  res <- rbindlist(lapply(unique(mk$chr), function(ch) {
    s <- mk[chr == ch]; if (nrow(s) < 10) return(NULL)
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    W <- matrix(NA_integer_, length(starts), ncol(G))
    for (k in seq_along(starts)) {
      rr <- which(s$ppos >= starts[k] & s$ppos < starts[k] + WIN * 1e6)
      if (!length(rr)) next
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    co <- apply(W, 2, function(v) { keep <- !is.na(v); if (sum(keep) < 2) return(NA_integer_)
      kk <- which(keep); sum(v[kk[-1]] != v[kk[-length(kk)]]) })
    data.table(dataset = f, chr = ch, co = mean(co, na.rm = TRUE))
  }))
  out[[f]] <- res
  cat("done", f, "\n")
}
R <- rbindlist(out)
R[, family := sub("_WATDE.*$", "", dataset)]
R <- R[, .(co = round(mean(co, na.rm = TRUE), 1)), by = .(family, chr)]
fwrite(R, file.path(TAB, "Table10_CO_by_chr_family.csv"))
cat("rows", nrow(R), "\n")
