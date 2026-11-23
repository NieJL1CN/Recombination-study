## 83_ld_by_CO_level.R -- LD decay and blocks for high vs low CO individuals (within family)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, whole := R1 + R2a + C + R2b + R3]
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
r2pair <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) < 15) return(NA_real_); cor(x[ok], y[ok])^2 }
res <- list(); blk <- list()
for (ds in COH) {
  G <- read_geno(ds); if (is.null(G)) next
  fam <- sub("_WATDE.*$", "", ds)
  m <- MK[dataset == ds]; m <- m[match(intersect(rownames(G), m$marker), marker)]
  if (nrow(m) < 100) next
  m <- m[order(match(chr, CHR), cM)]
  G <- G[m$marker, , drop = FALSE]
  y <- CO[dataset == ds][match(colnames(G), id)]$whole
  if (sum(!is.na(y)) < 40) next
  q <- quantile(y, c(0.33, 0.67), na.rm = TRUE)
  grp <- ifelse(y <= q[1], "low", ifelse(y >= q[2], "high", NA))
  for (g in c("high", "low")) {
    idx <- which(grp == g); if (length(idx) < 8) next
    for (ch in CHR) {
      sel <- which(m$chr == ch); if (length(sel) < 10) next
      sub <- G[sel, idx, drop = FALSE]; pos <- m$pos[sel]
      set.seed(1); np <- 400
      a <- sample(seq_along(sel), np, TRUE); b <- sample(seq_along(sel), np, TRUE)
      k <- a != b; a <- a[k]; b <- b[k]
      dd <- abs(pos[a] - pos[b])
      rr <- mapply(function(u, v) r2pair(sub[u, ], sub[v, ]), a, b)
      for (bin in list(c(0,5), c(5,10), c(10,20), c(20,50))) {
        s <- !is.na(rr) & dd >= bin[1] & dd < bin[2]
        if (sum(s) < 20) next
        res[[length(res) + 1]] <- data.table(family = fam, group = g, bin = paste0(bin[1], "-", bin[2]),
          mean_r2 = mean(rr[s]), n = sum(s))
      }
    }
    for (li in idx) {
      gg <- G[, li]; ok <- !is.na(gg); if (sum(ok) < 100) next
      runs <- rle(gg[ok])
      blk[[length(blk) + 1]] <- data.table(family = fam, group = g, line = colnames(G)[li],
        n_runs = length(runs$lengths), mean_run = mean(runs$lengths))
    }
  }
}
R <- rbindlist(res); B <- rbindlist(blk)
fwrite(R, file.path(TAB, "Table71_LD_by_CO_level.csv"))
fwrite(B, file.path(TAB, "Table72_blocks_by_CO_level.csv"))
print(R[, .(mean_r2 = round(mean(mean_r2), 3), n = .N), by = .(group, bin)][order(bin, group)])
print(B[, .(n_runs = round(mean(n_runs), 1), mean_run = round(mean(mean_run), 1), n = .N), by = group])
W <- dcast(R[, .(mean_r2 = mean(mean_r2)), by = .(family, bin, group)], family + bin ~ group, value.var = "mean_r2")
W <- W[!is.na(high) & !is.na(low)]
cat("\npaired Wilcoxon LD (high vs low): p =", signif(wilcox.test(W$high, W$low, paired = TRUE)$p.value, 3), " n =", nrow(W), "\n")
Bb <- B[, .(n_runs = mean(n_runs), mean_run = mean(mean_run)), by = .(family, group)]
Bw <- dcast(Bb, family ~ group, value.var = "n_runs"); Bw <- Bw[!is.na(high) & !is.na(low)]
cat("block number p =", signif(wilcox.test(Bw$high, Bw$low, paired = TRUE)$p.value, 3), " n =", nrow(Bw), "\n")
