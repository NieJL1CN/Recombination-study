## 00b_clean_tail.R -- drop trailing tie clusters that inflate a chromosome end
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived")
for (f in c("maps_min1.rds", "maps.rds")) {
  L <- readRDS(file.path(DER, f))
  CM <- as.data.table(L$consensus)
  CM[, mx := max(cM), by = .(panel, chr)]
  CM[, nmax := sum(cM == mx), by = .(panel, chr)]
  CM[, prev := max(cM[cM < mx], na.rm = TRUE), by = .(panel, chr)]
  bad <- unique(CM[nmax >= 5 & (mx - prev) > 20, .(panel, chr, nmax, mx, prev)])
  cat("\n", f, ": ", nrow(bad), " groups\n", sep = "")
  if (nrow(bad)) {
    print(bad)
    kill <- CM[bad, on = .(panel, chr)][cM == mx]
    cat("  removing", nrow(kill), "markers\n")
    CM <- CM[!marker %in% kill$marker]
    CM[, c("mx", "nmax", "prev") := NULL]
    L$consensus <- CM
    L$markers <- L$markers[marker %in% CM$marker | is.na(L$markers$g_rank)]
    saveRDS(L, file.path(DER, f))
  }
}