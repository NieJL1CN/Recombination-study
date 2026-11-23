suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
P <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P[, family := sub("_WATDE.*$", "", dataset)]
D <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
res <- list()
for (i in seq_len(nrow(D))) {
  d <- P[family == D$family[i] & trait == D$trait[i] & chr == D$chr[i]]
  if (!nrow(d)) next
  k <- which.max(d$lod)
  res[[i]] <- data.table(family = D$family[i], trait = D$trait[i], chr = D$chr[i],
    pos = d$pos[k], lod = d$lod[k], r2 = d$r2[k], eff = d$eff[k], thr = d$thr[k],
    direction = D$direction[i])
}
R <- rbindlist(res)
setorder(R, -r2)
R[, rank_r2 := seq_len(.N)]
R[, r2_pct := round(100 * r2, 1)]
print(R[, .(rank_r2, family, trait, chr, lod = round(lod,2), r2_pct, eff = round(eff,1), direction)])
fwrite(R, file.path(TAB, "Table94_QTL_effect_ranking.csv"))
