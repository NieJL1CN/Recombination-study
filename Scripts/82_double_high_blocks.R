## 82_double_high_blocks.R -- block fragmentation vs number of high-recombination alleles
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
norm <- function(x) gsub("[.-][0-9][ABD]$", "", x)
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
res <- list()
for (ds in COH) {
  G <- read_geno(ds); if (is.null(G)) next
  fam <- sub("_WATDE.*$", "", ds)
  m <- MK[dataset == ds]
  m <- m[match(intersect(rownames(G), m$marker), marker)]
  if (nrow(m) < 100) next
  m <- m[order(match(chr, CHR), cM)]
  G <- G[m$marker, , drop = FALSE]
  ## high-allele count per line across significant QTL available in this family
  hi <- rep(0L, ncol(G)); navail <- 0L
  for (i in seq_len(nrow(DIR))) {
    dd <- MK[family == fam & chr == DIR$chr[i]]
    dd <- dd[marker %in% rownames(G)]
    if (!nrow(dd)) next
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    g <- G[pk, ]
    highA <- DIR$direction[i] != "Watkins"
    hi <- hi + as.integer(if (highA) g == 0 else g == 1)
    navail <- navail + 1L
  }
  ## block statistics per line
  for (li in seq_len(ncol(G))) {
    g <- G[, li]; ok <- !is.na(g)
    if (sum(ok) < 100) next
    runs <- rle(g[ok])
    res[[length(res) + 1]] <- data.table(family = fam, dataset = ds, line = colnames(G)[li],
      n_high = hi[li], n_qtl_available = navail,
      n_runs = length(runs$lengths), mean_run = mean(runs$lengths))
  }
  cat(fam, "QTL available:", navail, " lines:", ncol(G), "\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table70_double_high_blocks.csv"))
cat("\nper-family summary (mean n_runs by n_high):\n")
print(R[, .(n_runs = round(mean(n_runs), 1), mean_run = round(mean(mean_run), 1), n = .N),
        by = .(family, n_high)][order(family, n_high)][n_high <= 3])
cat("\npooled by n_high (>=2 vs 0):\n")
print(R[, .(n_runs = round(mean(n_runs), 1), mean_run = round(mean(mean_run), 1), n = .N), by = n_high][order(n_high)])
A <- R[n_high == 0]$n_runs; B <- R[n_high >= 2]$n_runs
if (length(A) > 5 && length(B) > 5) cat("Wilcoxon n_runs (>=2 vs 0): p =", signif(wilcox.test(B, A)$p.value, 3), "\n")
