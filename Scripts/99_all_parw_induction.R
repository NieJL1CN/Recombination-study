## 99_all_parw_induction.R -- CS orthologs + meiosis induction for all ParW QTL genes
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); WD <- file.path(ROOT, "Source/WebData/blastdb")
H <- fread(file.path(TAB, "Table95_all_candidates_wheat.tsv"), header = FALSE,
           col.names = c("gene_id","target","pident","alen","evalue","bits"))
H[, acc := sub("^[a-z]+\\|([^|]+)\\|.*$", "\\1", target)]
X <- fread(file.path(WD, "wheat_uniprot_xref.tsv"), sep = "\t", quote = "")
setnames(X, c("acc","gene_names","protein_name","cs"))
X <- X[cs != ""]; X[, cs_id := sub(";.*$", "", cs)][, cs_id := sub("\\..*$", "", cs_id)]
setorder(H, gene_id, -bits)
T1 <- merge(H, X[, .(acc, cs_id)], by = "acc", all.x = TRUE)
T1 <- T1[!is.na(cs_id)][, head(.SD, 1), by = gene_id]
G <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
G <- G[grepl("^ParW", family)]
M <- merge(G, T1[, .(gene_id, cs_id, pident)], by = "gene_id", all.x = TRUE)
fwrite(M, file.path(TAB, "Table95_parw_genes_CS.csv"))
cat("ParW genes:", nrow(M), " mapped:", sum(!is.na(M$cs_id)), " families:", uniqueN(M$family), "\n")
ids <- unique(M[!is.na(cs_id)]$cs_id)
cat("unique CS ids:", length(ids), "\n")
writeLines(ids, file.path(TAB, "Table95_cs_ids.txt"))
