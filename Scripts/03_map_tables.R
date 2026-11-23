## ============================================================================
## 03_map_tables.R -- per-map marker tables and per-chromosome statistics
## ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER  <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR  <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps.rds")); M <- L$markers; CM <- L$consensus

## ---- Table 1a: per map, one row per marker --------------------------------
## columns: marker, consensus position, physical position, then each family's
## original genetic distance
fam_wide <- dcast(M[!is.na(g_rank), .(dataset, marker, g_rank)],
                  marker ~ dataset, value.var = "g_rank")
for (P in unique(CM$panel)) {
  D <- CM[panel == P, .(chr, marker, cM)]
  D <- merge(D, unique(M[, .(marker, ppos)]), by = "marker")
  fams <- intersect(names(fam_wide), unique(M[panel == P]$dataset))
  D <- merge(D, fam_wide[, c("marker", fams), with = FALSE], by = "marker", all.x = TRUE)
  setorder(D, chr, cM)
  fwrite(D, file.path(TAB, paste0("Table1a_markers_", P, ".csv")))
}

## ---- Table 1b: per map x chromosome statistics ----------------------------
stat <- function(d) {
  d <- d[order(cM)]
  g <- diff(d$cM); gp <- diff(d$ppos) / 1e6
  data.table(n_marker = nrow(d),
             first_marker = d$marker[1], last_marker = d$marker[nrow(d)],
             phys_start_Mb = round(min(d$ppos) / 1e6, 2),
             phys_end_Mb = round(max(d$ppos) / 1e6, 2),
             phys_len_Mb = round((max(d$ppos) - min(d$ppos)) / 1e6, 2),
             gen_len_cM = round(max(d$cM) - min(d$cM), 1),
             density_per_Mb = round(nrow(d) / max((max(d$ppos) - min(d$ppos)) / 1e6, 0.01), 2),
             density_per_cM = round(nrow(d) / max(max(d$cM) - min(d$cM), 0.01), 2),
             median_gap_cM = round(median(g), 2),
             max_gap_cM = round(max(g), 2),
             cM_per_Mb = round((max(d$cM) - min(d$cM)) / max((max(d$ppos) - min(d$ppos)) / 1e6, 0.01), 2))
}
S <- list()
for (P in unique(CM$panel)) {
  D <- merge(CM[panel == P, .(chr, marker, cM)], unique(M[, .(marker, ppos)]), by = "marker")
  S[[P]] <- D[, c(list(panel = P, level = "chromosome", chr = chr), stat(.SD)), by = chr][, chr := NULL]
  per_chr <- D[, c(list(panel = P, level = "chromosome", chr = chr), stat(.SD)), by = chr][, chr := NULL]
  S[[paste0(P, "_genome")]] <- cbind(panel = P, level = "genome", chr = "all",
    per_chr[, .(n_marker = sum(n_marker),
                phys_len_Mb = round(sum(phys_len_Mb), 2),
                gen_len_cM = round(sum(gen_len_cM), 1),
                density_per_Mb = round(sum(n_marker) / sum(phys_len_Mb), 2),
                density_per_cM = round(sum(n_marker) / sum(gen_len_cM), 2),
                median_gap_cM = round(median(D[, diff(sort(cM)), by = chr]$V1), 2),
                max_gap_cM = round(max(D[, diff(sort(cM)), by = chr]$V1), 2),
                cM_per_Mb = round(sum(gen_len_cM) / sum(phys_len_Mb), 2))])
  per_chr[, grp := sub("^[1-7]", "", chr)]
  S[[paste0(P, "_group")]] <- per_chr[, .(panel = P, level = "genome_group", group = grp,
                n_marker = sum(n_marker), phys_len_Mb = round(sum(phys_len_Mb), 2),
                gen_len_cM = round(sum(gen_len_cM), 1),
                density_per_Mb = round(sum(n_marker) / sum(phys_len_Mb), 2),
                cM_per_Mb = round(sum(gen_len_cM) / sum(phys_len_Mb), 2)), by = grp][, grp := NULL]
}
ST <- rbindlist(S, fill = TRUE)
setcolorder(ST, c("panel", "level", "chr", "group"))
fwrite(ST, file.path(TAB, "Table1b_map_statistics.csv"))
print(ST[level == "genome"])
cat("tables done\n")