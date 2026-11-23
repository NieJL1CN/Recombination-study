## 20_co_corrected.R -- CO with >=2 markers per window and short-segment filtering
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

run_co <- function(f, WIN, STEP, MINMARK = 2, MINGAP = 10) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  tot <- numeric(ncol(G))
  for (ch in unique(mk$chr)) {
    s <- mk[chr == ch]; if (nrow(s) < 10) next
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    W <- matrix(NA_integer_, length(starts), ncol(G)); MID <- starts + WIN * 1e6 / 2
    for (k in seq_along(starts)) {
      rr <- which(s$ppos >= starts[k] & s$ppos < starts[k] + WIN * 1e6)
      if (length(rr) < MINMARK) next
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    for (j in seq_len(ncol(G))) {
      v <- W[, j]; id <- which(!is.na(v)); if (length(id) < 2) next
      st <- v[id]
      ## run-length encoding; merge runs shorter than MINGAP (Mb) into the flanking state
      r <- rle(st)
      endbp <- c(MID[id][cumsum(r$lengths)])
      startbp <- c(MID[id][1], endbp[-length(endbp)])
      len <- (endbp - startbp) / 1e6
      keep <- r$lengths
      short <- which(len < MINGAP)
      if (length(short)) {
        for (k in short) {
          if (k == 1 || k == length(r$values)) next
          if (r$values[k - 1] == r$values[k + 1]) r$values[k] <- r$values[k - 1]
        }
      }
      r2 <- rle(r$values)
      tot[j] <- tot[j] + (length(r2$values) - 1)
    }
  }
  data.table(dataset = f, co = tot)
}
res <- list()
for (f in COH) { r <- run_co(f, 20, 10); if (!is.null(r)) { r[, win := "20/10"]; res[[f]] <- r } }
R <- rbindlist(res)
R[, family := sub("_WATDE.*$", "", dataset)]
S <- R[, .(mean_co = round(mean(co), 1)), by = .(dataset, family)][order(-mean_co)]
fwrite(S, file.path(TAB, "Table14_CO_corrected.csv"))
print(S[1:6]); cat("overall mean:", round(mean(S$mean_co), 1), "\n")
