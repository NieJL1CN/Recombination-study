## 68_mlh1_polymorphism.R -- MLH1-region polymorphism per family
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes"); DER <- file.path(ROOT, "Source/Derived")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, ppos)]
W <- MK[chr == "3D" & abs(ppos/1e6 - 616.22) <= 3]
cat("markers near MLH1:", nrow(W), " families:", uniqueN(W$dataset), "\n")
res <- list()
for (ds in unique(W$dataset)) {
  f <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv"))
  if (!file.exists(f)) next
  x <- fread(f, header = FALSE)
  hdr <- as.character(unlist(x[1, ]))
  dat <- x[-(1:3)]
  setnames(dat, hdr)
  setnames(dat, 1, "sample")
  mk <- hdr[-1]
  keep <- mk[mk %in% W[dataset == ds]$marker]
  sel <- dat[, c("sample", keep), with = FALSE]
  g <- as.matrix(sel[, -1])
  codes <- unique(as.vector(g))
  nA <- colSums(g == "A", na.rm = TRUE); nB <- colSums(g == "B", na.rm = TRUE)
  res[[ds]] <- data.table(dataset = ds,
    n_markers = ncol(g),
    n_poly = sum(nA > 3 & nB > 3),
    frac_missing = round(mean(g %in% c("H", "N", "", NA)), 3))
}
R <- rbindlist(res)
R[, family := sub("_WATDE.*$", "", dataset)]
R[, watde := sub("^[^_]+_([^_]+)_.*$", "\\1", dataset)]
setorder(R, -n_poly)
print(R[, .(family, watde, n_markers, n_poly, frac_missing)])