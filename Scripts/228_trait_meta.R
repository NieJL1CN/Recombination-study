## 228_trait_meta.R -- cluster per-season QTL into trait meta-QTL
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
Q <- fread(file.path(TAB, "Table175_trait_QTL_by_season.csv"))
setorder(Q, trait, chr, pos)
Q[, grp := cumsum(c(1, diff(pos) > 20 | chr[-1] != chr[-.N])), by = trait]
M <- Q[, .(chr = chr[1], pos = weighted.mean(pos, pmax(lod, 0.1)), n_qtl = .N,
           n_families = uniqueN(family), n_seasons = uniqueN(season),
           max_lod = round(max(lod), 2), start = min(pos), end = max(pos)), by = .(trait, grp)]
M[, grp := NULL]
fwrite(M, file.path(TAB, "Table176_trait_metaQTL.csv"))
print(M[, .N, by = trait]); print(head(M, 10))
