## 54_paragon_candidates.R -- candidate genes in CO-QTL intervals (Paragon annotation)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.63.chr.gff3.gz"), "r")
lines <- readLines(con, warn = FALSE)
close(con)
lines <- lines[grepl("\tgene\t", lines)]
sp <- tstrsplit(lines, "\t")
gff <- data.table(chr = sp[[1]], start = as.numeric(sp[[4]]), end = as.numeric(sp[[5]]),
                  strand = sp[[7]], attr = sp[[9]])
gff[, gene_id := sub(".*gene_id=([^;]+).*", "\\1", attr)]
gff[, biotype := sub(".*biotype=([^;]+).*", "\\1", attr)]
gff[, attr := NULL]
cat("genes:", nrow(gff), " chroms:", uniqueN(gff$chr), "\n")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
PK <- fread(file.path(TAB, "Table30b_CO_QTL_suggestive_peaks.csv"))
PK[, phys := NA_real_]
for (i in seq_len(nrow(PK))) {
  d <- MK[family == PK$family[i] & chr == PK$chr[i]]
  if (nrow(d)) PK$phys[i] <- d[which.min(abs(cM - PK$pos[i])), ppos]
}
PK <- PK[!is.na(phys)]
WIN <- 5e6
out <- list()
for (i in seq_len(nrow(PK))) {
  ch <- PK$chr[i]; s0 <- PK$phys[i] - WIN; e0 <- PK$phys[i] + WIN
  g <- gff[chr == ch & end >= s0 & start <= e0]
  if (!nrow(g)) next
  out[[i]] <- data.table(
    qtl_id = paste0(PK$trait[i], "-", PK$family[i], "-", ch, "-", round(PK$pos[i], 1)),
    family = PK$family[i], trait = PK$trait[i], chr = ch,
    peak_Mb = round(PK$phys[i] / 1e6, 2),
    interval_Mb = paste0(round(s0 / 1e6, 1), "-", round(e0 / 1e6, 1)),
    gene_id = g$gene_id, start = g$start, end = g$end, strand = g$strand, biotype = g$biotype)
}
G <- rbindlist(out)
fwrite(G, file.path(TAB, "Table45_CO_QTL_candidate_genes.csv"))
SUM <- G[, .(n_genes = uniqueN(gene_id)), by = .(qtl_id, family, trait, chr, peak_Mb, interval_Mb)]
fwrite(SUM, file.path(TAB, "Table45b_CO_QTL_gene_counts.csv"))
cat("candidate rows:", nrow(G), " QTL:", uniqueN(G$qtl_id), " genes/QTL median:", median(SUM$n_genes), "\n")
print(SUM[order(-n_genes)][1:8])