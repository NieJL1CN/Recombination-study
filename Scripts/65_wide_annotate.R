## 65_wide_annotate.R -- annotate the 4,176 genes in CO-QTL +-5 Mb windows
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
H <- fread(file.path(TAB, "Table53b_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, protein := sub("^sp\\|[^|]+\\|[^ ]+ ", "", desc)]
H[, species := sub(".*OS=([^ ]+) OX=.*", "\\1", protein)]
H[, protein := sub(" OS=.*$", "", protein)]
setorder(H, gene_id, -bits)
TOPH <- H[, head(.SD, 1), by = gene_id]
KW <- "RAD51|DMC1|SPO11|MSH[0-9]|MLH[0-9]|HEI10|ZYP1|ASY[0-9]|REC8|SCC[0-9]|FANCM|FIGL1|TOP2|CENH3|SMC[0-9]|MRE11|RAD50|NBS1|COM1|SIX1|synaptonemal|chiasma|crossover|meiotic|meiosis"
TOPH[, recomb := grepl(KW, protein, ignore.case = TRUE)]
G <- fread(file.path(TAB, "Table45_CO_QTL_candidate_genes.csv"))
A <- merge(G, TOPH[, .(gene_id, protein, species, evalue, recomb)], by = "gene_id", all.x = TRUE)
fwrite(A, file.path(TAB, "Table54_wide_annotation.csv"))
cat("rows:", nrow(A), " unique genes:", uniqueN(A$gene_id),
    " annotated:", uniqueN(A[!is.na(protein)]$gene_id),
    " recomb-flagged genes:", uniqueN(A[recomb == TRUE]$gene_id), "\n")
R <- unique(A[recomb == TRUE, .(gene_id, chr, protein, species, evalue)])
print(R[order(-evalue)][1:25])
