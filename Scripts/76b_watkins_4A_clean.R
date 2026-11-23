## 76b_watkins_4A_clean.R -- clean rebuild of the Watkins 4A candidate table
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); DOC <- file.path(ROOT, "Docs")
G <- fread(file.path(TAB, "Table56_significant_QTL_genes_5Mb.csv"))
S <- G[family == "ParW729" & chr == "4A"]
H <- fread(file.path(TAB, "Table56c_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, protein := sub(" OS=.*$", "", sub("^[^|]+\\|[^|]+\\|[^ ]+ ", "", desc))]
H[, species := sub(".*OS=([^ ]+) OX=.*", "\\1", desc)]
H[, entry := tstrsplit(target, "|", fixed = TRUE)[[3]]]
H[, sym := sub("_.*$", "", entry)]
H[, code := sub(".*_", "", entry)]
spmap <- c(ARATH="Arabidopsis thaliana", ORYSJ="rice (japonica)", ORYSI="rice (indica)",
           HORVU="barley", BOVIN="cow", DANRE="zebrafish", LOTJA="Lotus japonicus", YEAST="yeast")
H[, species2 := ifelse(code %in% names(spmap), spmap[code], code)]
setorder(H, gene_id, -bits)
H1 <- H[, head(.SD, 1), by = gene_id]
S <- merge(S, H1[, .(gene_id, protein, species2, evalue, pident)], by = "gene_id", all.x = TRUE)
S[, gene_Mb := round(gene_start/1e6, 2)]
S[, dist := round((gene_start - peak_Mb*1e6)/1e6, 2)]
setorder(S, gene_start)
out <- S[, .(gene_id, chr, gene_Mb, dist, family, trait, direction, lod, protein, species2, evalue, pident)]
setnames(out, c("gene_id","chr","gene_Mb","dist_Mb","family","trait","direction","lod",
                "swissprot_protein","species","evalue","identity_pct"))
fwrite(out, file.path(TAB, "Table64_Watkins_QTL_4A_candidates.csv"))
cat("genes:", nrow(out), " with hits:", sum(!is.na(out$swissprot_protein)), "\n")
cat("peak Mb:", unique(S$peak_Mb), " LOD:", unique(S$lod), "\n\n")
for (i in seq_len(nrow(out))) cat(sprintf("%-26s %7.2f (%+.2f)  %-45s %s\n",
  out$gene_id[i], out$gene_Mb[i], out$dist_Mb[i],
  ifelse(is.na(out$swissprot_protein[i]), "no hit", substr(out$swissprot_protein[i], 1, 45)),
  ifelse(is.na(out$species[i]), "", out$species[i])))
