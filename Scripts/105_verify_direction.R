## 105_verify_direction.R -- verify QTL direction by direct group means at peak marker
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]; CO[, middle := R2a + R2b]; CO[, centromere := C]; CO[, whole := R1 + R2a + C + R2b + R3]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
out <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]; tr <- DIR$trait[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  best <- 0; res <- NULL
  for (cand in ds_all) {
    fp <- file.path(GEN, paste0(cand, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
    x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
    mk <- colnames(x)[-1]
    dd <- MK[dataset == cand & chr == ch]; dd <- dd[marker %in% mk]
    if (nrow(dd) < 10) next
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    gt <- g[[pk]]
    y <- CO[match(g[[1]], id)][[tr]]
    ok <- gt %in% c("A","B") & !is.na(y)
    if (sum(ok) < 30) next
    mA <- mean(y[ok & gt == "A"]); mB <- mean(y[ok & gt == "B"])
    res <- data.table(family = fam, chr = ch, trait = tr, marker = pk,
      nA = sum(ok & gt == "A"), nB = sum(ok & gt == "B"),
      mean_A = round(mA, 2), mean_B = round(mB, 2), diff_B_minus_A = round(mB - mA, 2),
      eff_table = round(DIR$eff[i], 2), direction_table = DIR$direction[i])
    break
  }
  if (!is.null(res)) out[[length(out)+1]] <- res
}
R <- rbindlist(out)
print(R)
fwrite(R, file.path(TAB, "Table100_direction_verification.csv"))
cat("\nconsistent:", sum(sign(R$diff_B_minus_A) == sign(R$eff_table)), "of", nrow(R), "\n")
