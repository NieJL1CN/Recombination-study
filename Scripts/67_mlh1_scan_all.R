## 67_mlh1_scan_all.R -- LOD at the MLH1 locus for every family + polymorphism check
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
MLH_MB <- 616.22
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
## consensus markers within +-3 Mb of MLH1 on 3D
W <- MK[chr == "3D" & abs(ppos/1e6 - MLH_MB) <= 3]
cat("markers within 3 Mb:", nrow(W), " cM range:", paste(round(range(W$cM), 1), collapse = "-"), "\n")
cm_lo <- min(W$cM); cm_hi <- max(W$cM)
## LOD per family at that interval
P <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P[, family := sub("_WATDE.*$", "", dataset)]
S <- P[chr == "3D" & pos >= cm_lo - 1 & pos <= cm_hi + 1,
       .(lod = max(lod), thr = thr[which.max(lod)], eff = eff[which.max(lod)], pos = pos[which.max(lod)]),
       by = .(family, trait)]
S[, direction := ifelse(eff > 0, "Watkins", "Paragon")]
S[, watde := MK[family == family][1, sub("^[^_]+_([^_]+)_.*$", "\\1", dataset)]]
setorder(S, -lod)
fwrite(S, file.path(TAB, "Table55_MLH1_locus_LOD_by_family.csv"))
cat("\nTop families at MLH1 locus:\n")
print(head(S[, .(family, watde, trait, pos_cM = round(pos,1), lod = round(lod,2),
                 thr = round(thr,2), direction)], 15))
cat("\nfamilies with LOD > thr:", S[lod > thr, .N], " rows\n")
print(S[lod > thr, .(family, trait, lod = round(lod,2), thr = round(thr,2), direction)])
