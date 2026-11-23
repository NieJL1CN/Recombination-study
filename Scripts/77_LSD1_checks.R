## 77_LSD1_checks.R -- subfamily, map compression, developmental association
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
## (1) map compression: markers near 7D:393.91 Mb
W <- MK[chr == "7D" & abs(ppos/1e6 - 393.91) <= 2]
cat("7D markers within +-2 Mb of 393.91:\n")
print(unique(W[, .(marker, cM = round(cM, 2), ppos_Mb = round(ppos/1e6, 2))])[order(cM)])
cat("\ncM range in that window:", paste(round(range(W$cM), 2), collapse = "-"), "\n")
cat("7D map length:", round(max(MK[chr == "7D"]$cM), 1), "cM\n")
## (2) published HD/Ant QTL near 7D:393.91
P <- fread(file.path(TAB, "Table35_published_QTL_cohort.csv"))
P[, chr2 := sub("_.*$", "", chr)]
H <- P[trait %in% c("Ant_dto_day","Hd_dto_day") & chr2 == "7D" & !is.na(bp_peak)]
cat("\npublished HD QTL on 7D within 20 Mb of 393.91:\n")
print(H[abs(bp_peak/1e6 - 393.91) <= 20, .(population, trait, bp_Mb = round(bp_peak/1e6, 1),
      LOD = round(LOD, 2), diff_Mb = round(bp_peak/1e6 - 393.91, 1))])
## (3) peak marker allele vs heading date in ParW352
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))[family == "ParW352" & chr == "7D"]
cat("\nParW352 7D peak cM:", PK$pos, " LOD:", round(PK$lod, 2), "\n")
ds <- unique(MK[family == "ParW352" & chr == "7D"]$dataset)
d <- MK[family == "ParW352" & chr == "7D"]
near <- d[which.min(abs(cM - PK$pos))]
cat("nearest marker:", near$marker, " cM:", round(near$cM, 2), " Mb:", round(near$ppos/1e6, 2), " dataset:", near$dataset, "\n")
