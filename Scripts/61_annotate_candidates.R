## 61_annotate_candidates.R -- parse DIAMOND hits and flag recombination genes
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
H <- fread(file.path(TAB, "Table48b_diamond_sprot.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, desc := sub("^sp\\|[^|]+\\|[^ ]+ ", "", desc)]
H[, desc := sub(" OS=.*$", "", desc)]
H[, gene_sym := sub(" .*$", "", desc)]
H[, species := sub(".*OS=([^=]+) OX=.*", "\\1", target)]
TOPH <- H[order(gene_id, -bits)][, head(.SD, 1), by = gene_id]
KW <- "recombin|meiotic|meiosis|crossover|RAD51|DMC1|SPO11|MSH[0-9]|MLH[0-9]|HEI10|ZYP1|ASY[0-9]|REC8|SCC[0-9]|FANCM|FIGL1|TOP2|CENH3|SMC[0-9]|MRE11|RAD50|NBS1|COM1|SIX1|synaptonemal|chiasma"
TOPH[, recomb := grepl(KW, desc, ignore.case = TRUE)]
G <- fread(file.path(TAB, "Table47_significant_QTL_candidate_genes_1Mb.csv"))
A <- merge(TOPH, G, by = "gene_id", all.x = TRUE)
fwrite(A, file.path(TAB, "Table50_candidate_gene_annotation.csv"))
cat("annotated:", nrow(A), " recombination-flagged:", sum(A$recomb), "\n")
print(A[order(-recomb, -bits), .(gene_id, chr, peak_Mb, desc, species, evalue, recomb)][1:25])
