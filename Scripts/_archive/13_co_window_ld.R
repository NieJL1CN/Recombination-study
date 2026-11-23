## ============================================================================
## 13_co_window_ld.R -- CO counting with the LD-derived window (5.5 Mb, step 2.75)
##   marker-based CO and window-based CO for every individual
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- 5.5; STEP <- 2.75
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)

M <- readRDS(file.path(DER, "maps_min1.rds"))$markers
M <- M[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

fs <- unique(M$dataset)
out <- list()
for (f in fs) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv"))
  if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- data.table(marker = rownames(G))
  mk <- merge(mk, pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- merge(mk, B, by = "chr")
  mk[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
        fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
  mk <- mk[order(match(chr, CHR), ppos)]

  ## marker-based switches
  ids <- colnames(G)
  mark <- rbindlist(lapply(seq_along(ids), function(j) {
    v <- G[, j]
    data.table(dataset = f, id = ids[j],
               co_marker = sum(diff(v) != 0, na.rm = TRUE))
  }))

  ## window-based switches
  win <- rbindlist(lapply(unique(mk$chr), function(ch) {
    s <- mk[chr == ch]
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    ws <- data.table(x = starts); ws[, xend := x + WIN * 1e6]
    ws[, reg := s$region[findInterval(x + WIN * 1e6 / 2, s$ppos) + 1]]
    ws <- ws[!is.na(reg)]
    if (nrow(ws) < 2) return(NULL)
    W <- matrix(NA_integer_, nrow(ws), ncol(G))
    for (k in seq_len(nrow(ws))) {
      rr <- which(s$ppos >= ws$x[k] & s$ppos < ws$xend[k])
      if (!length(rr)) next
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    rbindlist(lapply(seq_along(ids), function(j) {
      v <- W[, j]; keep <- !is.na(v); if (sum(keep) < 2) return(NULL)
      kk <- which(keep)
      sw <- kk[-1][v[kk[-1]] != v[kk[-length(kk)]]]
      if (!length(sw)) return(NULL)
      data.table(dataset = f, id = ids[j], chr = ch, region = ws$reg[sw], co_win = 1)
    }))
  }))
  Wtot <- win[, .(co_win = .N), by = .(dataset, id)]
  A <- merge(mark, Wtot, by = c("dataset", "id"), all.x = TRUE)
  A[is.na(co_win), co_win := 0L]
  out[[f]] <- A
  cat("done", f, "\n")
}
CO <- rbindlist(out, fill = TRUE)
fwrite(CO, file.path(TAB, "Table9_CO_by_individual.csv"))
cat("individuals:", nrow(CO), "\n")
print(CO[, .(co_marker = round(mean(co_marker),1), co_win = round(mean(co_win),1)), by = dataset][order(-co_win)][1:5])