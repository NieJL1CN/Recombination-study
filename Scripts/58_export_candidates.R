## 58_export_candidates.R -- export protein sequences for candidate genes
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
G <- fread(file.path(TAB, "Table47_significant_QTL_candidate_genes_1Mb.csv"))
M <- fread(file.path(TAB, "Table30_CO_metaQTL.csv"))[n_families >= 2]
SEL <- G[chr %in% M$chr]
keep <- unique(SEL$gene_id)
cat("candidate genes in meta-QTL chromosomes:", length(keep), "\n")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.pep.all.fa.gz"), "r")
res <- list(); cur <- NULL; buf <- character()
repeat {
  ln <- readLines(con, n = 50000, warn = FALSE)
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
cat("sequences found:", length(res), "of", length(keep), "\n")
fa <- file.path(TAB, "Table48_candidate_proteins.faa")
con2 <- file(fa, "w")
for (g in names(res)) cat(">", g, "\n", res[[g]], "\n", sep = "", file = con2)
close(con2)
cat("written", fa, "\n")
