## 66_mlh1_families.R -- which families carry a CO-QTL at the MLH1 locus (3D)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
## MLH1 gene coordinates from Paragon GFF
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.63.chr.gff3.gz"), "r")
lines <- readLines(con, warn = FALSE); close(con)
lines <- lines[grepl("\\tgene\\t", lines)]
sp <- tstrsplit(lines, "\\t")
g <- data.table(chr = sp[[1]], start = as.numeric(sp[[4]]), end = as.numeric(sp[[5]]), attr = sp[[9]])
g[, gene_id := sub(".*gene_id=([^;]+).*", "\\1", attr)]
MLH <- g[gene_id %in% c("TraesPARA_EIv1.0_1171450","TraesPARA_EIv1.0_1171550","TraesPARA_EIv1.0_1171500")]
cat("MLH1 genes:\n"); print(MLH[, .(gene_id, chr, start, end, Mb = round(start/1e6, 2))])
mlh_mb <- mean(MLH$start) / 1e6
cat("MLH1 centre (Mb):", round(mlh_mb, 2), "\n\n")
## all CO-QTL peaks with physical position
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
MK[, watde := sub("^[^_]+_([^_]+)_.*$", "\\1", dataset)]
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
PK[, phys := NA_real_]
for (i in seq_len(nrow(PK))) {
  d <- MK[family == PK$family[i] & chr == PK$chr[i]]
  if (nrow(d)) PK$phys[i] <- d[which.min(abs(cM - PK$pos[i])), ppos]
}
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))
PK <- merge(PK, DIR[, .(family, trait, chr, pos, direction, significant)], by = c("family","trait","chr","pos"), all.x = TRUE)
SEL <- PK[chr == "3D" & !is.na(phys) & abs(phys/1e6 - mlh_mb) <= 15]
cat("CO-QTL within 15 Mb of MLH1:\n")
print(SEL[, .(family, watde = MK[family == family][1, watde], trait, pos_cM = round(pos,1),
              phys_Mb = round(phys/1e6,2), dist_Mb = round(phys/1e6 - mlh_mb, 2), lod = round(lod,2),
              direction, significant)])
## also check phenotype QTL (Table2c) near MLH1
PH <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))[map == "consensus"]
PH[, phys := NA_real_]
for (i in seq_len(nrow(PH))) {
  d <- MK[family == PH$family[i] & chr == PH$chr[i]]
  if (nrow(d)) PH$phys[i] <- d[which.min(abs(cM - PH$pos[i])), ppos]
}
P2 <- PH[chr == "3D" & !is.na(phys) & abs(phys/1e6 - mlh_mb) <= 15]
cat("\nphenotype QTL within 15 Mb of MLH1:\n")
print(P2[, .(family, trait = tr_short, pos = round(pos,1), phys_Mb = round(phys/1e6,2), lod = round(lod,2))])
