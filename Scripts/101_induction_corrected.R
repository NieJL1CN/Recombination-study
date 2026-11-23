## 101_induction_corrected.R -- log2(TPM+1) induction (correct scale)
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
parse1 <- function(pattern) {
  fs <- list.files(TAB, pattern = pattern, full.names = TRUE); out <- list()
  for (f in fs) { j <- try(fromJSON(f, simplifyVector = FALSE), silent = TRUE)
    if (inherits(j, "try-error") || is.null(j$data$results)) next
    for (g in j$data$results) for (p in g$points)
      out[[length(out)+1]] <- data.table(cs_id = g$gene_id, stage = p$label, value = as.numeric(p$value)) }
  rbindlist(out)
}
mei <- parse1("^expr_all_PRJEB5029"); tis <- parse1("^expr_all_PRJEB5314")
e4  <- fread(file.path(TAB, "expr_4B_genes.json")); t4 <- fread(file.path(TAB, "expr_4B_tissue.json"))
cat("scales are TPM; applying log2(x+1)\n")
ind <- function(mei_dt, veg_dt) {
  M <- mei_dt[, .(mei = max(log2(value + 1), na.rm = TRUE)), by = cs_id]
  V <- veg_dt[grepl("leaf|root|stem", stage, ignore.case = TRUE),
              .(veg = max(log2(value + 1), na.rm = TRUE)), by = cs_id]
  merge(M, V, by = "cs_id", all = TRUE)[, induction := mei - veg][]
}
## all ParW
D <- ind(mei, tis)
G <- fread(file.path(TAB, "Table95_parw_genes_CS.csv"))
A <- unique(fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))[, .(gene_id, swiss = protein, recomb)])
R <- merge(merge(G, D, by = "cs_id", all.x = TRUE), A, by = "gene_id", all.x = TRUE)
fwrite(R, file.path(TAB, "Table96_all_parw_induction.csv"))
cat("genes:", nrow(R), " with expr:", sum(!is.na(R$induction)), "\n")
cat("\n== top 12 induced per QTL (log2) ==\n")
print(R[!is.na(induction)][order(family, chr, -induction)][, head(.SD, 3), by = .(family, chr)][
  , .(family, chr, gene_id, ind = round(induction,1), swiss = substr(ifelse(is.na(swiss), "NA", swiss), 1, 42))])
cat("\n== recombination-keyword genes ==\n")
print(R[recomb == TRUE & !is.na(induction)][order(-induction)][
  , .(family, chr, gene_id, ind = round(induction,1), swiss = substr(swiss, 1, 45))])
## 4B subset
D4 <- ind(e4, t4)
G4 <- fread(file.path(TAB, "Table89_4B_QTL_genes_CS.csv"))
A4 <- merge(merge(G4, D4, by = "cs_id", all.x = TRUE), A, by = "gene_id", all.x = TRUE)
fwrite(A4, file.path(TAB, "Table93_4B_genes_induction.csv"))
cat("\n== 4B QTL genes (corrected) ==\n")
print(A4[!is.na(induction)][order(-induction)][1:8, .(gene_id, cs_id, ind = round(induction,2),
      swiss = substr(ifelse(is.na(swiss), "NA", swiss), 1, 40))])
