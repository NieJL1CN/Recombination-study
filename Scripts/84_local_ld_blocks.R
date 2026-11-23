## 84_local_ld_blocks.R -- LD and ancestry blocks on the QTL chromosome only
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
r2pair <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) < 15) return(NA_real_); cor(x[ok], y[ok])^2 }
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
LD <- list(); BL <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  got <- FALSE
  for (cand in ds_all) {
    G <- read_geno(cand); if (is.null(G)) next
    dd <- MK[dataset == cand & chr == ch]; dd <- dd[marker %in% rownames(G)]
    if (nrow(dd) < 20) next
    dd <- dd[order(cM)]
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    gt <- G[pk, ]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == 0 else gt == 1
    non <- if (highA) gt == 1 else gt == 0
    if (sum(car, na.rm = TRUE) < 10 || sum(non, na.rm = TRUE) < 10) next
    for (g in c("carrier", "non")) {
      idx <- if (g == "carrier") which(car) else which(non)
      sub <- G[dd$marker, idx, drop = FALSE]
      pos <- dd$ppos / 1e6
      set.seed(1); np <- 400
      a <- sample(seq_along(dd$marker), np, TRUE); b <- sample(seq_along(dd$marker), np, TRUE)
      k <- a != b; a <- a[k]; b <- b[k]
      dist <- abs(pos[a] - pos[b])
      rr <- mapply(function(u, v) r2pair(sub[u, ], sub[v, ]), a, b)
      for (bin in list(c(0,2), c(2,5), c(5,10), c(10,20), c(20,50))) {
        s <- !is.na(rr) & dist >= bin[1] & dist < bin[2]
        if (sum(s) < 15) next
        LD[[length(LD) + 1]] <- data.table(family = fam, qtl_chr = ch, group = g,
          bin = paste0(bin[1], "-", bin[2]), mean_r2 = mean(rr[s]), n = sum(s))
      }
      for (li in idx) {
        gg <- G[dd$marker, li]; ok <- !is.na(gg)
        if (sum(ok) < 20) next
        runs <- rle(gg[ok])
        BL[[length(BL) + 1]] <- data.table(family = fam, qtl_chr = ch, group = g,
          line = colnames(G)[li], n_runs = length(runs$lengths),
          mean_run_markers = mean(runs$lengths),
          mean_run_Mb = mean(runs$lengths) * mean(diff(pos[ok])[diff(pos[ok]) >= 0], na.rm = TRUE))
      }
    }
    got <- TRUE
    break
  }
  cat(fam, ch, ifelse(got, "ok", "skipped"), "\n")
}
L2 <- rbindlist(LD); B2 <- rbindlist(BL)
fwrite(L2, file.path(TAB, "Table73_local_LD_by_genotype.csv"))
fwrite(B2, file.path(TAB, "Table74_local_blocks_by_genotype.csv"))
cat("\nlocal LD (QTL chromosome):\n")
print(L2[, .(mean_r2 = round(mean(mean_r2), 3), n = .N), by = .(group, bin)][order(bin, group)])
cat("\nlocal blocks:\n")
print(B2[, .(n_runs = round(mean(n_runs), 1), run_markers = round(mean(mean_run_markers), 1), n = .N), by = group])
W <- dcast(L2[, .(mean_r2 = mean(mean_r2)), by = .(family, bin, group)], family + bin ~ group, value.var = "mean_r2")
W <- W[!is.na(carrier) & !is.na(non)]
if (nrow(W) > 3) {
  cat("\npaired Wilcoxon local LD: p =", signif(wilcox.test(W$carrier, W$non, paired = TRUE)$p.value, 3), " n =", nrow(W), " mean diff =", round(mean(W$carrier - W$non), 4), " negatives =", sum(W$carrier < W$non), "\n")
}
Bb <- B2[, .(n_runs = mean(n_runs)), by = .(family, group)]
Bw <- dcast(Bb, family ~ group, value.var = "n_runs"); Bw <- Bw[!is.na(carrier) & !is.na(non)]
if (nrow(Bw) > 3) cat("paired Wilcoxon local block number: p =", signif(wilcox.test(Bw$carrier, Bw$non, paired = TRUE)$p.value, 3), " n =", nrow(Bw), "\n")
