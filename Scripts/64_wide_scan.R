## 64_wide_scan.R -- export all genes within +-5 Mb of CO-QTL and annotate
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
G <- fread(file.path(TAB, "Table45_CO_QTL_candidate_genes.csv"))
keep <- unique(G$gene_id)
cat("genes to export:", length(keep), "\n")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.pep.all.fa.gz"), "r")
res <- list(); cur <- NULL; buf <- character()
repeat {
  ln <- readLines(con, n = 100000, warn = FALSE)
  if (!length(ln)) break
  for (x in ln) {
    if (startsWith(x, ">")) {
      if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
      gid <- sub(".*gene:([^ ]+).*", "\\1", x)
      cur <- if (gid %in% keep) gid else NULL
      buf <- character()
    } else if (!is.null(cur)) buf <- c(buf, trimws(x))
  }
}
if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
close(con)
cat("sequences:", length(res), "\n")
con2 <- file(file.path(TAB, "Table53_wide_candidate_proteins.faa"), "w")
for (g in names(res)) cat(">", g, "\n", res[[g]], "\n", sep = "", file = con2)
close(con2)
