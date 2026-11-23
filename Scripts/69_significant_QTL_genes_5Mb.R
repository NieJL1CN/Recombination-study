## 69_significant_QTL_genes_5Mb.R -- candidate genes in SIGNIFICANT CO-QTL (+-5 Mb)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
## significant CO-QTL only
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))
SIG <- DIR[significant == TRUE, .(family, trait, chr, pos, lod, eff, direction)]
cat("significant QTL:", nrow(SIG), "\n")
MK <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank),
        .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
SIG[, phys := NA_real_]
for (i in seq_len(nrow(SIG))) {
  d <- MK[family == SIG$family[i] & chr == SIG$chr[i]]
  if (nrow(d)) SIG$phys[i] <- d[which.min(abs(cM - SIG$pos[i])), ppos]
}
SIG <- SIG[!is.na(phys)]
cat("with physical position:", nrow(SIG), "\n")
## genes
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.63.chr.gff3.gz"), "r")
lines <- readLines(con, warn = FALSE); close(con)
lines <- lines[grepl("\\tgene\\t", lines)]
sp <- tstrsplit(lines, "\\t")
gff <- data.table(chr = sp[[1]], start = as.numeric(sp[[4]]), end = as.numeric(sp[[5]]), attr = sp[[9]])
gff[, gene_id := sub(".*gene_id=([^;]+).*", "\\1", attr)]
gff[, attr := NULL]
WIN <- 5e6
out <- list()
for (i in seq_len(nrow(SIG))) {
  ch <- SIG$chr[i]; s0 <- SIG$phys[i] - WIN; e0 <- SIG$phys[i] + WIN
  g <- gff[chr == ch & end >= s0 & start <= e0]
  if (!nrow(g)) next
  out[[i]] <- data.table(qtl_id = paste0(SIG$trait[i], "-", SIG$family[i], "-", ch),
    family = SIG$family[i], trait = SIG$trait[i], chr = ch, peak_Mb = round(SIG$phys[i]/1e6, 2),
    lod = round(SIG$lod[i], 2), direction = SIG$direction[i],
    gene_id = g$gene_id, gene_start = g$start, gene_end = g$end,
    dist_Mb = round((g$start - SIG$phys[i]) / 1e6, 2))
}
G <- rbindlist(out)
fwrite(G, file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
cat("gene rows:", nrow(G), " unique genes:", uniqueN(G$gene_id), " QTL:", uniqueN(G$qtl_id), "\n")
print(G[, .(n_genes = uniqueN(gene_id)), by = qtl_id][order(-n_genes)])
## export proteins
keep <- unique(G$gene_id)
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
cat("proteins exported:", length(res), "\n")
con2 <- file(file.path(TAB, "Table56b_proteins.faa"), "w")
for (g in names(res)) cat(">", g, "\n", res[[g]], "\n", sep = "", file = con2)
close(con2)
