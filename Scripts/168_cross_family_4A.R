## 168_cross_family_4A.R -- test the Paragon 4A allele across families
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds")); M <- L$markers
PH <- fread(file.path(ROOT, "Tables/Table22c_CO_rate_individual.csv"))
out <- list()
for (fam in unique(sub("_WATDE.*$", "", M$dataset))) {
  cand <- unique(M[grepl(paste0("^", fam, "_"), dataset), dataset])
  if (!length(cand)) next
  ds <- cand[which.max(sapply(cand, function(z) nrow(M[dataset == z & !is.na(ppos)])))]
  d <- M[dataset == ds & chr == "4A" & !is.na(ppos)]
  if (!nrow(d)) next
  best <- d[which.min(abs(ppos - 530.9e6))]
  f <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv"))
  if (!file.exists(f)) next
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  ids <- x[[1]][-(1:2)]
  g <- x[-(1:2), -1, drop = FALSE]
  if (!best$marker %in% colnames(g)) next
  gt <- trimws(g[[best$marker]])
  gt[gt %in% c("A","a","1","AA")] <- "A"; gt[gt %in% c("B","b","2","BB")] <- "B"
  w <- PH[dataset == ds]
  w[, id := sub(paste0(fam, "_"), "", id)]
  id2 <- sub(paste0(fam, "_"), "", ids)
  win <- w[match(id2, w$id)]
  co <- rowSums(win[, .(R1, R3)], na.rm = TRUE)
  ok <- gt %in% c("A","B") & !is.na(co)
  if (sum(ok) < 40) next
  a <- co[ok & gt == "A"]; b <- co[ok & gt == "B"]
  tt <- t.test(b, a)
  out[[length(out)+1]] <- data.table(family = fam, marker = best$marker, Mb = round(best$ppos/1e6,1),
    nA = length(a), nB = length(b), meanA = round(mean(a),2), meanB = round(mean(b),2),
    effect = round(mean(b) - mean(a),2), p = tt$p.value)
}
R <- rbindlist(out)[order(p)]
R[, direction := ifelse(effect > 0, "B higher", "A higher")]
print(R)
cat("\nfamilies with B (Paragon-like) higher:", sum(R$effect > 0), "of", nrow(R), "\n")
cat("families with P<0.05:", sum(R$p < 0.05), "\n")
fwrite(R, file.path(ROOT, "Tables/Table151_cross_family_4A.csv"))
