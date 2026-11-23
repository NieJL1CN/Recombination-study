## 95_4B_genes_expression.R -- all genes in the 4B QTL: CS orthologs + meiosis expression
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); WD <- file.path(ROOT, "Source/WebData/blastdb")
G <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))[family == "ParW579" & chr == "4B"]
cat("genes in 4B QTL window:", nrow(G), "\n")
H <- fread(file.path(TAB, "Table56c_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, acc := sub("^[a-z]+\\|([^|]+)\\|.*$", "\\1", target)]
X <- fread(file.path(WD, "wheat_uniprot_xref.tsv"), sep = "\t", quote = "")
setnames(X, c("acc","gene_names","protein_name","cs"))
X <- X[cs != ""]; X[, cs_id := sub(";.*$", "", cs)][, cs_id := sub("\\..*$", "", cs_id)]
X[, chr := sub("^TraesCS([0-9][ABD]).*$", "\\1", cs_id)]
R <- merge(H, X[, .(acc, cs_id, chr)], by = "acc", all.x = TRUE)
setorder(R, gene_id, -bits)
T1 <- R[, head(.SD, 1), by = gene_id]
M <- merge(G[, .(gene_id, chr, peak_Mb, direction)], T1[, .(gene_id, cs_id, pident, evalue)],
           by = "gene_id", all.x = TRUE)
fwrite(M, file.path(TAB, "Table89_4B_QTL_genes_CS.csv"))
cat("with CS ortholog:", sum(!is.na(M$cs_id)), "of", nrow(M), "\n")
print(M)
