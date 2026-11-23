## 74_homeolog_structure.R -- gene models via transcript-parent mapping
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.63.chr.gff3.gz"), "r")
lines <- readLines(con, warn = FALSE); close(con)
lines <- lines[!startsWith(lines, "#")]
sp <- tstrsplit(lines, "\t")
g <- data.table(chr = sp[[1]], type = sp[[3]], start = as.numeric(sp[[4]]), end = as.numeric(sp[[5]]),
                strand = sp[[7]], attr = sp[[9]])
GENES <- c("TraesPARA_EIv1.0_1336720","TraesPARA_EIv1.0_1262790","TraesPARA_EIv1.0_1437800",
           "TraesPARA_EIv1.0_2580590","TraesPARA_EIv1.0_2301160","TraesPARA_EIv1.0_2423550")
## gene rows
gg <- g[type == "gene" & grepl(paste(GENES, collapse = "|"), attr)]
gg[, gid := sub(".*gene:([^;]+).*", "\\1", attr)]
## mRNA rows -> gene
mm <- g[type == "mRNA"]
mm[, tx := sub(".*ID=([^;]+);.*", "\\1", attr)]
mm[, gid := sub(".*Parent=gene:([^;]+).*", "\\1", attr)]
mm <- mm[gid %in% GENES]
tx2gene <- mm[, .(tx, gid)]
## exons / CDS
ex <- g[type %in% c("exon","CDS")]
ex[, tx := sub(".*Parent=([^;]+).*", "\\1", attr)]
ex <- merge(ex, tx2gene, by = "tx")
S <- ex[, .(n_exon = sum(type == "exon"), n_cds = sum(type == "CDS"),
            cds_len = sum((end - start + 1)[type == "CDS"])), by = gid]
S <- merge(S, gg[, .(gid, chr, start, end, strand)], by = "gid", all.x = TRUE)
setorder(S, gid)
fwrite(S, file.path(TAB, "Table61_homeolog_gene_models.csv"))
print(S)
