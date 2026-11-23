## 169_bootstrap_qtl.R
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds")); M <- L$markers
PK <- fread(file.path(ROOT, "Tables/Table20b_CO_QTL_peaks.csv"))
PK[, family := sub("_WATDE.*$", "", dataset)]
PH <- fread(file.path(ROOT, "Tables/Table22c_CO_rate_individual.csv"))
set.seed(20261211)
out <- list()
for (i in seq_len(nrow(PK))) {
  fam <- PK$family[i]; ch <- PK$chr[i]; ds <- PK$dataset[i]
  d <- M[dataset == ds & chr == ch & !is.na(ppos)]
  if (!nrow(d)) next
  best <- d[which.min(abs(g_rank - PK$pos[i]))]
  f <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv"))
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  ids <- x[[1]][-(1:2)]
  g <- x[-(1:2), -1, drop = FALSE]
  if (!length(best$marker)) next
  if (!best$marker %in% colnames(g)) next
  gt <- trimws(g[[best$marker]])
  gt[gt %in% c("A","a","1","AA")] <- 0; gt[gt %in% c("B","b","2","BB")] <- 1
  gt <- as.numeric(gt)
  w <- PH[dataset == ds]; w[, id := sub(paste0(fam, "_"), "", id)]
  id2 <- sub(paste0(fam, "_"), "", ids)
  wi <- w[match(id2, w$id)]
  trait <- PK$trait[i]
  co <- if (trait == "distal") rowSums(wi[, .(R1, R3)], na.rm = TRUE) else
        if (trait == "centromere") wi$C else
        rowSums(wi[, .(R1, R2a, C, R2b, R3)], na.rm = TRUE)
  ok <- !is.na(gt) & !is.na(co)
  y <- co[ok]; xg <- gt[ok]; n <- length(y)
  if (n < 40) next
  eff <- mean(y[xg == 1]) - mean(y[xg == 0])
  lod <- function(y, x) { fit <- lm(y ~ x); s <- summary(fit); -log10(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE))/2 }
  L0 <- lod(y, xg)
  bs <- replicate(1000, { idx <- sample.int(n, n, TRUE); lod(y[idx], xg[idx]) })
  be <- replicate(1000, { idx <- sample.int(n, n, TRUE); mean(y[idx][xg[idx] == 1]) - mean(y[idx][xg[idx] == 0]) })
  out[[i]] <- data.table(qtl = paste0(fam, "-", ch), trait = trait, n = n, effect = round(eff,2), LOD = round(L0,2),
    boot_LOD_median = round(median(bs),2), boot_LOD_2.5 = round(quantile(bs,0.025),2), boot_LOD_97.5 = round(quantile(bs,0.975),2),
    sign_consistency = round(mean(sign(be) == sign(eff)),3), prop_LOD_gt2 = round(mean(bs > 2),3))
}
R <- rbindlist(out)
print(R)
fwrite(R, file.path(ROOT, "Tables/Table152_bootstrap_qtl.csv"))
