## 91_map_to_CS.R -- map Paragon candidates to Chinese Spring (IWGSC) gene IDs
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); WD <- file.path(ROOT, "Source/WebData/blastdb")
H <- fread(file.path(TAB, "Table83b_CS_hits.tsv"), header = FALSE,
           col.names = c("paragon", "target", "pident", "alen", "evalue", "bits"))
H[, acc := sub("^[a-z]+\\|([^|]+)\\|.*$", "\\1", target)]
X <- fread(file.path(WD, "wheat_uniprot_xref.tsv"), sep = "\t", quote = "")
setnames(X, c("acc","gene_names","protein_name","cs"))
X <- X[cs != ""]
X[, cs_id := sub(";.*$", "", cs)]
X[, cs_id := sub("\\..*$", "", cs_id)]
X[, chr := sub("^TraesCS([0-9][ABD]).*$", "\\1", cs_id)]
R <- merge(H, X[, .(acc, cs_id, chr, gene_names, protein_name)], by = "acc", all.x = TRUE)
setorder(R, paragon, -bits)
R[, rank := seq_len(.N), by = paragon]
fwrite(R, file.path(TAB, "Table84_Paragon_to_CS_mapping.csv"))
top <- R[rank <= 3]
for (p in unique(top$paragon)) {
  cat("\n==", p, "==\n")
  print(top[paragon == p, .(rank, cs_id, chr, pident, evalue = signif(evalue, 2), gene_names)])
}
