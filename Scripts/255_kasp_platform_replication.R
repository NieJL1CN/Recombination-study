## 255_kasp_platform_replication.R
## Cross-platform replication of the crossover QTL with the KASP panel.
## The crossover phenotype is the AXIOM-derived one (Table22c); the KASP panel is
## used only to genotype the peak region, so the test is not circular.
## Output: Tables/Table184_KASP_platform_replication.csv (marker-level results),
##         summarised as Table S21 in the manuscript.
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
TAB <- file.path(ROOT, "Tables")
MAP <- as.data.table(readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))$markers)
MAP <- MAP[!is.na(ppos), .(marker, chr, ppos)][!duplicated(marker)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]
CO[, middle := R2a + R2b]
CO[, centromere := C]
CO[, whole := R1 + R2a + C + R2b + R3]
qtl <- data.table(
  family = c("ParW579","ParW729","ParW729","ParW352","ParW579","ParW049","ParW007","ParW145","RobW145"),
  chr    = c("3B","1B","4A","7D","4B","7B","6A","5D","6B"),
  peak   = c(582.0, 330.0, 530.9, 393.6, 167.3, 3.7, 568.5, 544.6, 459.8),
  trait  = c("distal","distal","distal","distal","distal","centromere","whole","whole","whole"))
out <- list()
for (k in seq_len(nrow(qtl))) {
  fam <- qtl$family[k]; tchr <- qtl$chr[k]; mb <- qtl$peak[k]; tr <- qtl$trait[k]
  f <- list.files(GEN, pattern = paste0("^", fam, "_.*KASP_rqtl_genotypes\\.csv$"), full.names = TRUE)
  if (!length(f)) { cat(fam, tchr, ": no KASP file\n"); next }
  x <- fread(f[1], colClasses = "character"); ids <- x[[1]][-c(1, 2)]
  sub <- MAP[marker %in% names(x) & chr == tchr & abs(ppos/1e6 - mb) <= 60]
  if (!nrow(sub)) { cat(fam, tchr, ": no KASP marker within 60 Mb\n"); next }
  d <- CO[grepl(paste0("^", fam, "_"), id), c("id", tr), with = FALSE]
  setnames(d, tr, "pheno")
  best <- NULL
  for (i in seq_len(nrow(sub))) {
    m <- sub$marker[i]; gv <- x[[m]][-c(1, 2)]
    gv[gv %in% c("A","a","AA")] <- "0"; gv[gv %in% c("B","b","BB")] <- "1"
    gv[!gv %in% c("0","1")] <- NA
    dd <- data.table(id = ids, geno = as.numeric(gv))[!is.na(geno)]
    dd <- merge(dd, d, by = "id")
    a <- dd[geno == 0, pheno]; b <- dd[geno == 1, pheno]
    if (length(a) < 8 || length(b) < 8) next
    p <- tryCatch(t.test(b, a)$p.value, error = function(e) NA)
    r <- data.table(family = fam, chr = tchr, trait = tr, peak = mb, marker = m,
                    Mb = sub$ppos[i]/1e6, n0 = length(a), n1 = length(b),
                    mean0 = mean(a), mean1 = mean(b), effect = mean(b) - mean(a), p = p)
    if (is.null(best) || (!is.na(p) && p < best$p)) best <- r
    out[[length(out) + 1]] <- r
  }
  if (!is.null(best)) cat(sprintf("%-8s %-2s %-10s best %.1f Mb (%s): effect=%+.2f p=%.4f\n",
     fam, tchr, tr, best$Mb, best$marker, best$effect, best$p))
}
res <- rbindlist(out)
fwrite(res, file.path(TAB, "Table184_KASP_platform_replication.csv"))
cat("Table184 written:", nrow(res), "marker-level rows\n")
