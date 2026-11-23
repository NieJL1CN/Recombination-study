## 57_paragon_candidates_tight.R -- candidate genes (+-1 Mb) around significant CO-QTL
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
## genes from GFF
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.63.chr.gff3.gz"), "r")
lines <- readLines(con, warn = FALSE); close(con)
lines <- lines[grepl("\\tgene\\t", lines)]
sp <- tstrsplit(lines, "\\t")
gff <- data.table(chr = sp[[1]], start = as.numeric(sp[[4]]), end = as.numeric(sp[[5]]),
                  strand = sp[[7]], attr = sp[[9]])
gff[, gene_id := sub(".*gene_id=([^;]+).*", "\\1", attr)]
## significant CO-QTL with physical position
PK <- fread(file.path(TAB, "Table30b_CO_QTL_suggestive_peaks.csv"))
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))
SIG <- DIR[significant == TRUE, .(family, trait, chr, pos, lod, eff, direction)]
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
SIG[, phys := NA_real_]
for (i in seq_len(nrow(SIG))) {
  d <- MK[family == SIG$family[i] & chr == SIG$chr[i]]
  if (nrow(d)) SIG$phys[i] <- d[which.min(abs(cM - SIG$pos[i])), ppos]
}
SIG <- SIG[!is.na(phys)]
WIN <- 1e6
out <- list()
for (i in seq_len(nrow(SIG))) {
  ch <- SIG$chr[i]; s0 <- SIG$phys[i] - WIN; e0 <- SIG$phys[i] + WIN
  g <- gff[chr == ch & end >= s0 & start <= e0]
  if (!nrow(g)) next
  out[[i]] <- data.table(qtl_id = paste0(SIG$trait[i], "-", SIG$family[i], "-", ch),
    chr = ch, peak_Mb = round(SIG$phys[i] / 1e6, 2), direction = SIG$direction[i],
    gene_id = g$gene_id, start = g$start, end = g$end, strand = g$strand)
}
G <- rbindlist(out)
fwrite(G, file.path(TAB, "Table47_significant_QTL_candidate_genes_1Mb.csv"))
cat("QTL:", uniqueN(G$qtl_id), " gene rows:", nrow(G), " unique genes:", uniqueN(G$gene_id), "\n")
print(G[, .(n_genes = uniqueN(gene_id)), by = qtl_id][order(-n_genes)])
