## 100_all_parw_induction_summary.R -- induction for all ParW QTL genes
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
parse_expr <- function(pattern) {
  fs <- list.files(TAB, pattern = pattern, full.names = TRUE)
  out <- list()
  for (f in fs) {
    j <- try(fromJSON(f, simplifyVector = FALSE), silent = TRUE)
    if (inherits(j, "try-error") || is.null(j$data$results)) next
    for (g in j$data$results) for (p in g$points)
      out[[length(out)+1]] <- data.table(cs_id = g$gene_id, stage = p$label, value = as.numeric(p$value))
  }
  rbindlist(out)
}
mei <- parse_expr("^expr_all_PRJEB5029")
tis <- parse_expr("^expr_all_PRJEB5314")
cat("meiosis rows:", nrow(mei), " tissue rows:", nrow(tis), "\n")
M <- mei[, .(mei_max = max(value, na.rm = TRUE)), by = cs_id]
V <- tis[grepl("leaf|root|stem", stage, ignore.case = TRUE),
         .(veg_max = max(value, na.rm = TRUE)), by = cs_id]
D <- merge(unique(M), V, by = "cs_id", all = TRUE)
D[, induction := mei_max - veg_max]
G <- fread(file.path(TAB, "Table95_parw_genes_CS.csv"))
A <- fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))
A <- unique(A[, .(gene_id, swiss = protein, recomb)])
R <- merge(G, D, by = "cs_id", all.x = TRUE)
R <- merge(R, A, by = "gene_id", all.x = TRUE)
fwrite(R, file.path(TAB, "Table96_all_parw_induction.csv"))
cat("genes with expression:", sum(!is.na(R$induction)), "of", nrow(R), "\n\n")
cat("top 15 induced overall:\n")
print(R[!is.na(induction)][order(-induction)][1:15, .(family, chr, gene_id, cs_id,
      mei = round(mei_max,1), veg = round(veg_max,1), ind = round(induction,1), swiss)])
cat("\nrecombination-keyword genes with induction:\n")
print(R[recomb == TRUE & !is.na(induction)][order(-induction)][1:15, .(family, chr, gene_id,
      ind = round(induction,1), swiss)])
cat("\nper-QTL median induction and count induced >2:\n")
print(R[!is.na(induction), .(n = .N, med = round(median(induction),2), n_ind2 = sum(induction > 2)),
        by = .(family, chr)][order(-n_ind2)])
