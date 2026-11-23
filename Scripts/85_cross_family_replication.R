## 85_cross_family_replication.R -- test each QTL peak in all other families
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
P <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
norm <- function(x) gsub("[.-][0-9][ABD]$", "", x)
MK[, nm := norm(marker)]
out <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]; tr <- DIR$trait[i]
  d <- MK[family == fam & chr == ch]
  pk <- d[which.min(abs(cM - DIR$pos[i]))]$nm
  ## for every family: is this marker present, at what cM, LOD?
  for (f2 in unique(MK$family)) {
    dd <- MK[family == f2 & chr == ch & nm == pk]
    if (!nrow(dd)) { out[[length(out)+1]] <- data.table(qtl = paste0(tr,"-",fam,"-",ch), test_family = f2, tested = FALSE, lod = NA_real_); next }
    cm <- dd$cM[1]
    lod <- P[family == f2 & trait == tr & chr == ch & abs(pos - cm) <= 2]$lod
    out[[length(out)+1]] <- data.table(qtl = paste0(tr,"-",fam,"-",ch), test_family = f2,
      tested = TRUE, lod = if (length(lod)) max(lod) else NA_real_)
  }
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table75_cross_family_replication.csv"))
cat("QTLs:", uniqueN(R$qtl), " families:", uniqueN(R$test_family), "\n\n")
print(R[, .(share_tested = round(mean(tested), 2), max_LOD_other = round(max(lod, na.rm = TRUE), 2),
            n_sig_other = sum(lod > 3, na.rm = TRUE)), by = qtl][order(-max_LOD_other)])