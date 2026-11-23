suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
TAB <- file.path(ROOT, "Tables")
MAP <- as.data.table(readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))$markers)
MAP <- MAP[!is.na(ppos), .(marker, chr)][!duplicated(marker)]
files <- list.files(GEN, pattern = "_rqtl_genotypes.csv$", full.names = TRUE)
out <- list()
for (fp in files) {
  ds <- sub("_rqtl_genotypes.csv$", "", basename(fp))
  x <- fread(fp, colClasses = "character")
  G <- as.matrix(x[-c(1,2), -1, drop = FALSE])
  mk <- colnames(x)[-1]
  G[G %in% c("A","a","1","AA")] <- 0; G[G %in% c("B","b","2","BB")] <- 1
  G[!(G %in% c("0","1"))] <- NA
  storage.mode(G) <- "double"
  M <- MAP[match(mk, MAP$marker)]
  ok <- !is.na(M$chr)
  G <- G[, ok, drop = FALSE]; M <- M[ok]
  for (ch in unique(M$chr)) {
    sel <- which(M$chr == ch)
    n <- sum(rowSums(!is.na(G[, sel, drop = FALSE])) > 0)
    out[[length(out)+1]] <- data.table(dataset = ds, chr = ch, n_chr = n)
  }
}
NC <- rbindlist(out)
fwrite(NC, file.path(TAB, "Table182_n_RIL_per_chr.csv"))
E <- fread(file.path(TAB, "Table178_marker_intervals.csv"))
E <- merge(E, NC, by = c("dataset","chr"))
W <- 10e6
E[, b0 := floor(a/W)*W][, b1 := floor(b/W)*W]
B <- E[, .(bin = seq(b0[1], b1[1], by = W)), by = .(dataset, chr, a, b, d, n_chr)]
B[, w := pmax(0, pmin(bin + W, b) - pmax(bin, a)) / d / n_chr]
P <- B[, .(co = sum(w)), by = .(dataset, chr, bin)][, .(co = mean(co)), by = .(dataset, chr, bin)]
A <- P[, .(co = mean(co) * 100 / 10), by = .(chr, bin)]
fwrite(A, file.path(TAB, "Table179_CO_10Mb_marker.csv"))
cat("max", round(max(A$co*10),2), "\n")
print(A[chr %in% c("4D","5B") & bin/1e6 %in% c(10,50,100,150,200,250,300), .(chr, Mb = bin/1e6, co10 = round(co*10,2))])
