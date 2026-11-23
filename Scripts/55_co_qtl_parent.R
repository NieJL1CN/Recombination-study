## 55_co_qtl_parent.R -- which parent contributes the CO-increasing allele?
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
P <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P[, family := sub("_WATDE.*$", "", dataset)]
## peak per family x trait x chr
PK <- P[order(chr, pos), .SD[which.max(lod)], by = .(family, trait, chr)]
PK <- PK[, .(family, trait, chr, pos, lod, eff, thr, r2)]
PK[, significant := lod > thr]
## A allele (coded 0) is the non-Watkins parent: Paragon for ParW*, Robigus for RobW*
PK[, parentA := ifelse(grepl("^RobW", family), "Robigus", "Paragon")]
## empirical check: negative eff = B (Watkins) allele increases; positive = A allele
PK[, direction := ifelse(eff < 0, "Watkins", parentA)]
PK[, peak_Mb_note := NA_character_]
fwrite(PK, file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))
cat("QTL total:", nrow(PK), " significant:", sum(PK$significant), "\n\n")
cat("direction among significant QTL:\n")
print(PK[significant == TRUE, .N, by = .(trait, direction)][order(trait)])
cat("\ndirection among all suggestive QTL (lod >= 2):\n")
print(PK[lod >= 2, .N, by = direction])
cat("\noverall significant:", PK[significant == TRUE, .N, by = direction])