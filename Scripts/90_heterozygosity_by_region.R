## 90_heterozygosity_by_region.R -- residual heterozygosity by region and QTL genotype
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
M0 <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
M0[, family := sub("_WATDE.*$", "", dataset)]
M1 <- merge(M0, B, by = "chr")
M1[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
       fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
out <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(M1[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    fp <- file.path(GEN, paste0(cand, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
    x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
    lines <- g[[1]]; mk <- colnames(x)[-1]; M <- as.matrix(g[, -1, drop = FALSE])
    dd <- M1[dataset == cand]; dd <- dd[match(mk, dd$marker)]
    ok <- !is.na(dd$region)
    if (sum(ok) < 100) next
    pk <- dd[chr == ch][which.min(abs(cM - DIR$pos[i]))]$marker
    if (!length(pk)) next
    gt <- M[, match(pk, mk)]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == "A" else gt == "B"
    non <- if (highA) gt == "B" else gt == "A"
    if (sum(car) < 10 || sum(non) < 10) next
    Hm <- (M == "H")
    for (rg in c("R1", "R2a", "C", "R2b", "R3")) {
      sel <- which(ok & dd$region == rg)
      if (length(sel) < 20) next
      h <- rowSums(Hm[, sel, drop = FALSE])
      out[[length(out) + 1]] <- data.table(family = fam, region = rg, line = lines,
        group = ifelse(car, "carrier", ifelse(non, "non", NA)),
        frac_H = h / length(sel), n_marker = length(sel))
    }
    break
  }
}
R <- rbindlist(out); R <- R[!is.na(group)]
fwrite(R, file.path(TAB, "Table82_heterozygosity_by_region.csv"))
cat("H fraction by region (all):\n")
print(R[, .(frac_H = round(mean(frac_H), 4), n = .N), by = region][order(-frac_H)])
cat("\nH fraction by region x group:\n")
print(dcast(R[, .(frac_H = mean(frac_H)), by = .(region, group)], region ~ group, value.var = "frac_H"))
cat("\ncarrier vs non by region (paired across families):\n")
S <- R[, .(frac_H = mean(frac_H)), by = .(family, region, group)]
for (rg in c("R1", "R2a", "C", "R2b", "R3")) {
  W <- dcast(S[region == rg], family ~ group, value.var = "frac_H")
  W <- W[!is.na(carrier) & !is.na(non)]
  if (nrow(W) > 3) cat(rg, ": diff =", signif(mean(W$carrier - W$non), 4),
    " p =", signif(wilcox.test(W$carrier, W$non, paired = TRUE)$p.value, 3), " n =", nrow(W), "\n")
}
