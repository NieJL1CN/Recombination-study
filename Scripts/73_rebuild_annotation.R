## 73_rebuild_annotation.R -- rebuild Table57 with corrected parent labels
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
H <- fread(file.path(TAB, "Table56c_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, protein := sub(" OS=.*$", "", sub("^[^|]+\\|[^|]+\\|[^ ]+ ", "", desc))]
H[, species := sub(".*OS=([^ ]+) OX=.*", "\\1", desc)]
setorder(H, gene_id, -bits)
T1 <- H[, head(.SD, 1), by = gene_id]
G <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
A <- merge(G, T1[, .(gene_id, protein, species, evalue, bits)], by = "gene_id", all.x = TRUE)
KW <- "RAD51|DMC1|SPO11|MSH[0-9]|MLH[0-9]|HEI10|ZYP1|ASY[0-9]|REC8|SCC[0-9]|FANCM|FIGL1|TOP2|CENH3|SMC[0-9]|MRE11|RAD50|NBS1|COM1|SIX1|synaptonemal|chiasma|crossover|meiotic|meiosis|histone|helicase|ubiquitin"
A[, recomb := grepl(KW, protein, ignore.case = TRUE)]
fwrite(A, file.path(TAB, "Table57_significant_QTL_annotation.csv"))
cat("rows:", nrow(A), " flagged:", uniqueN(A[recomb == TRUE & !is.na(recomb)]$gene_id), "\n")
print(A[recomb == TRUE, .N, by = direction])
