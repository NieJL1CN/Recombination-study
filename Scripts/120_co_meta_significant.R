suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
## ONLY significant QTL (LOD > permutation threshold)
S <- fread(file.path(TAB, "Table20b_CO_QTL_peaks.csv"))
S[, family := sub("_WATDE.*$", "", dataset)]
setorder(S, chr, pos)
S[, cluster := {
  cl <- integer(.N); cur <- 1L
  for (i in seq_len(.N)) {
    if (i > 1 && pos[i] - pos[i-1] > 20) cur <- cur + 1L
    cl[i] <- cur
  }
  cl
}, by = chr]
M <- S[, .(n_families = uniqueN(family), n_traits = uniqueN(trait), n_qtl = .N,
           start = min(pos), end = max(pos), mean_lod = mean(lod), max_lod = max(lod),
           families = paste(sort(unique(family)), collapse = ";"),
           traits = paste(sort(unique(trait)), collapse = ";")),
       by = .(chr, cluster)]
M[, multi := n_families >= 2]
fwrite(M, file.path(TAB, "Table30_CO_metaQTL.csv"))
cat("significant QTL:", nrow(S), " clusters:", nrow(M), " multi-family:", sum(M$multi), "\n")
print(M[, .(chr, cluster, n_families, n_traits, n_qtl, start = round(start,1), end = round(end,1),
            max_lod = round(max_lod,2), multi)])
