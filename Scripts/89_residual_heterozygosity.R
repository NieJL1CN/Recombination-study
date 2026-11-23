## 89_residual_heterozygosity.R -- test whether high-CO carriers are fixed faster
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
read_raw <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  list(G = g[, -1, drop = FALSE], lines = g[[1]], markers = colnames(x)[-1])
}
out <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    r <- read_raw(cand); if (is.null(r)) next
    M <- as.matrix(r$G)
    dd <- MK[dataset == cand & chr == ch]
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    if (!pk %in% r$markers) next
    gt <- M[, match(pk, r$markers)]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == "A" else gt == "B"
    non <- if (highA) gt == "B" else gt == "A"
    if (sum(car) < 10 || sum(non) < 10) next
    nH <- rowSums(M == "H", na.rm = TRUE)
    nNA <- rowSums(matrix(M %in% c("N",""), nrow = nrow(M)), na.rm = TRUE)
    out[[length(out) + 1]] <- data.table(family = fam, chr = ch, line = r$lines,
      group = ifelse(car, "carrier", ifelse(non, "non", NA)),
      n_H = nH, frac_H = nH / ncol(M), n_missing = nNA)
    break
  }
}
R <- rbindlist(out); R <- R[!is.na(group)]
fwrite(R, file.path(TAB, "Table81_residual_heterozygosity.csv"))
S <- R[, .(frac_H = mean(frac_H), n = .N), by = .(family, group)]
print(dcast(S, family ~ group, value.var = "frac_H"))
W <- dcast(S, family ~ group, value.var = "frac_H"); W <- W[!is.na(carrier) & !is.na(non)]
cat("\npaired Wilcoxon frac_H (carrier vs non): p =", signif(wilcox.test(W$carrier, W$non, paired = TRUE)$p.value, 3),
    " mean diff =", signif(mean(W$carrier - W$non), 4), "\n")
cat("pairs where carrier has fewer H:", sum(W$carrier < W$non), "of", nrow(W), "\n")