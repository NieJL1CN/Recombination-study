## 97_induction_and_promoter.R -- meiosis induction + corrected promoter comparison
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
parse_expr <- function(f) {
  j <- fromJSON(f, simplifyVector = FALSE); out <- list()
  for (g in j$data$results) for (p in g$points)
    out[[length(out)+1]] <- data.table(cs_id = g$gene_id, stage = p$label, value = as.numeric(p$value))
  rbindlist(out)
}
mei <- parse_expr(file.path(TAB, "expr_4B_genes.json"))
tis <- parse_expr(file.path(TAB, "expr_4B_tissue.json"))
cat("meiosis stages:", paste(unique(mei$stage), collapse = ", "), "\n")
cat("tissue stages:", paste(unique(tis$stage), collapse = ", "), "\n\n")
veg <- tis[stage %in% c("leaf","root","stem","Leaf","Root","Stem") | grepl("leaf|root|stem", stage, ignore.case = TRUE)]
V <- veg[, .(veg_max = max(value, na.rm = TRUE)), by = cs_id]
M <- mei[, .(mei_max = max(value, na.rm = TRUE)), by = cs_id]
D <- merge(M, V, by = "cs_id", all = TRUE)
D[, induction := mei_max - veg_max]
G <- fread(file.path(TAB, "Table89_4B_QTL_genes_CS.csv"))
D <- merge(D, G[, .(cs_id, gene_id, gene_names)], by = "cs_id", all.x = TRUE)
setorder(D, -induction)
fwrite(D, file.path(TAB, "Table93_4B_genes_induction.csv"))
print(D[!is.na(induction)][1:12, .(gene_id, cs_id, mei_max = round(mei_max,1),
      veg_max = round(veg_max,1), induction = round(induction,2))])
## promoter comparison (both oriented 5'->3' of the gene)
rc <- function(s) {
  s <- toupper(gsub("[^ACGT]", "", s))
  chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
}
read_fa <- function(f) paste(readLines(f)[!startsWith(readLines(f), ">")], collapse = "")
par_p <- rc(read_fa(file.path(TAB, "Table91_Paragon_RECQL1_4B_promoter.fa")))
cs_p  <- rc(read_fa(file.path(TAB, "Table92_CS_RECQL1_4B_upstream_plus.fa")))
cat("\npromoter lengths: Paragon", nchar(par_p), " CS", nchar(cs_p), "\n")
km <- function(s, k = 8) { n <- nchar(s) - k + 1; v <- substring(s, 1:n, 1:n + k - 1); unique(v[!grepl("[^ACGT]", v)]) }
a <- km(par_p); b <- km(cs_p)
cat("8-mer Jaccard:", round(length(intersect(a,b))/length(union(a,b)), 3),
    " shared:", length(intersect(a,b)), "of", length(a), "\n")
gc <- function(s) sum(strsplit(s,"")[[1]] %in% c("G","C"))/nchar(s)
cat("GC Paragon:", round(gc(par_p),3), " CS:", round(gc(cs_p),3), "\n")
