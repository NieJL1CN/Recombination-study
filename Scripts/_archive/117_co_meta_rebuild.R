## 117_co_meta_rebuild.R -- meta-QTL = clustering of ALL QTL (singletons included)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
P <- fread(file.path(TAB, "Table30b_CO_QTL_suggestive_peaks.csv"))
setorder(P, chr, pos)
## cluster: merge peaks within 20 cM on the same chromosome (all families/traits)
P[, cluster := {
  cl <- integer(.N); cur <- 1L
  for (i in seq_len(.N)) {
    if (i > 1 && pos[i] - pos[i-1] > 20) cur <- cur + 1L
    cl[i] <- cur
  }
  cl
}, by = chr]
M <- P[, .(n_families = uniqueN(family), n_traits = uniqueN(trait), n_qtl = .N,
           start = min(pos), end = max(pos), mean_lod = mean(lod), max_lod = max(lod),
           families = paste(sort(unique(family)), collapse = ";"),
           traits = paste(sort(unique(trait)), collapse = ";")),
       by = .(chr, cluster)]
M[, multi := n_families >= 2]
fwrite(M, file.path(TAB, "Table30_CO_metaQTL.csv"))
cat("clusters total:", nrow(M), " with >=2 families:", sum(M$multi), "\n")
print(M[order(-n_families, -max_lod)][1:12, .(chr, cluster, n_families, n_traits, n_qtl,
       start = round(start,1), end = round(end,1), max_lod = round(max_lod,2), multi)])
cat("\nsingle-family clusters:", sum(!M$multi), "\n")