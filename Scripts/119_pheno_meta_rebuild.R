suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
Q <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))[map == "consensus"]
Q[, family := sub("_WATDE.*$", "", dataset)]
setorder(Q, chr, pos)
Q[, cluster := {
  cl <- integer(.N); cur <- 1L
  for (i in seq_len(.N)) {
    if (i > 1 && pos[i] - pos[i-1] > 20) cur <- cur + 1L
    cl[i] <- cur
  }
  cl
}, by = chr]
M <- Q[, .(n_families = uniqueN(family), n_traits = uniqueN(tr_short), n_qtl = .N,
           start = min(pos), end = max(pos), mean_lod = mean(lod), max_lod = max(lod),
           families = paste(sort(unique(family)), collapse = ";"),
           traits = paste(sort(unique(tr_short)), collapse = ";")),
       by = .(chr, cluster)]
M[, multi := n_families >= 2]
fwrite(M, file.path(TAB, "Table2f_metaQTL.csv"))
cat("phenotype meta clusters:", nrow(M), " multi:", sum(M$multi), " single:", sum(!M$multi), "\n")
print(M[order(-n_families, -max_lod)][1:10, .(chr, cluster, n_families, n_traits, n_qtl,
      start = round(start,1), end = round(end,1), max_lod = round(max_lod,2), multi)])
