suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
TAB <- file.path(ROOT, "Tables")
MAP <- as.data.table(readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))$markers)
MAP <- MAP[!is.na(ppos), .(marker, chr, ppos)][!duplicated(marker)]
files <- list.files(GEN, pattern = "_rqtl_genotypes.csv$", full.names = TRUE)
res <- list()
for (fp in files) {
  ds <- sub("_rqtl_genotypes.csv$", "", basename(fp))
  x <- fread(fp, colClasses = "character")
  G <- as.matrix(x[-c(1,2), -1, drop = FALSE])
  mk <- colnames(x)[-1]
  G[G %in% c("A","a","1","AA")] <- 0; G[G %in% c("B","b","2","BB")] <- 1
  G[!(G %in% c("0","1"))] <- NA
  storage.mode(G) <- "double"
  colnames(G) <- mk
  M <- MAP[match(mk, MAP$marker)]
  keep <- !is.na(M$chr) & !is.na(M$ppos)
  G <- G[, keep, drop = FALSE]; M <- M[keep]
  out <- list()
  for (i in seq_len(nrow(G))) {
    v <- G[i, ]
    for (ch in unique(M$chr)) {
      sel <- which(M$chr == ch); if (length(sel) < 3) next
      vv <- v[sel]; pp <- M$ppos[sel]
      ok <- which(!is.na(vv)); if (length(ok) < 3) next
      vv <- vv[ok]; pp <- pp[ok]
      chg <- which(diff(vv) != 0); if (!length(chg)) next
      a <- pp[chg]; b <- pp[chg + 1]; d <- b - a
      good <- d > 0; if (!any(good)) next
      out[[length(out)+1]] <- data.table(dataset = ds, chr = ch, a = a[good], b = b[good], d = d[good])
    }
  }
  res[[length(res)+1]] <- rbindlist(out)
  cat(basename(fp), "\n")
}
E <- rbindlist(res)
fwrite(E, file.path(TAB, "Table178_marker_intervals.csv"))
W <- 10e6
E[, b0 := floor(a/W)*W][, b1 := floor(b/W)*W]
B <- E[, .(bin = seq(b0[1], b1[1], by = W)), by = .(dataset, chr, a, b, d)]
B[, w := pmax(0, pmin(bin + W, b) - pmax(bin, a)) / d]
P <- B[, .(co = sum(w)), by = .(dataset, chr, bin)][, .(co = mean(co)), by = .(dataset, chr, bin)]
A <- P[, .(co = mean(co) * 100), by = .(chr, bin)]
fwrite(A, file.path(TAB, "Table179_CO_10Mb_marker.csv"))
cat("intervals", nrow(E), "bins", nrow(A), "max", round(max(A$co),1), "\n")
print(A[chr == "5B" & bin/1e6 >= 220 & bin/1e6 <= 290])
