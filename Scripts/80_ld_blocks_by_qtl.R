## 80_ld_blocks_by_qtl.R -- LD decay and ancestry-block fragmentation by QTL genotype
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
CM <- L$consensus[, .(marker, chr, panel, pos = cM)][!duplicated(marker)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv"))
  if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  G
}
r2pair <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); if (sum(ok) < 20) return(NA_real_)
  cor(x[ok], y[ok])^2
}
res <- list(); blk <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  if (!length(ds_all)) next
  ds <- NA_character_; G <- NULL
  for (cand in ds_all) {
    gg <- read_geno(cand)
    if (is.null(gg)) next
    dd <- MK[dataset == cand & chr == ch]
    dd <- dd[marker %in% rownames(gg)]
    if (!nrow(dd)) next
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    if (length(pk)) { ds <- cand; G <- gg; break }
  }
  if (is.null(G)) next
  C <- CM[marker %in% rownames(G)]
  pnl <- C[, .N, by = panel][order(-N)]$panel[1]; C <- C[panel == pnl]
  mk <- C[order(match(chr, CHR), pos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  ## peak marker
  d <- MK[dataset == ds & chr == ch]
  peak <- d[which.min(abs(cM - DIR$pos[i]))]$marker
  if (!peak %in% rownames(G)) next
  gt <- G[peak, ]
  highA <- DIR$direction[i] != "Watkins"
  car <- if (highA) gt == 0 else gt == 1
  non <- if (highA) gt == 1 else gt == 0
  cat(fam, ch, "peak", peak, " carriers:", sum(car, na.rm = TRUE), " non:", sum(non, na.rm = TRUE), "\n")
  ## LD by group, per chromosome distance bins
  for (grp in c("carrier", "non")) {
    idx <- if (grp == "carrier") which(car) else which(non)
    if (length(idx) < 8) next
    for (ch2 in CHR) {
      sel <- which(mk$chr == ch2)
      if (length(sel) < 10) next
      sub <- G[sel, idx, drop = FALSE]
      pos <- mk$pos[sel]
      np <- 300
      set.seed(1)
      a <- sample(seq_along(sel), np, TRUE); b <- sample(seq_along(sel), np, TRUE)
      keep <- a != b
      a <- a[keep]; b <- b[keep]
      dd <- abs(pos[a] - pos[b])
      rr <- mapply(function(u, v) r2pair(sub[u, ], sub[v, ]), a, b)
      if (!any(!is.na(rr))) next
      for (bin in list(c(0,5), c(5,10), c(10,20), c(20,50))) {
        s <- !is.na(rr) & dd >= bin[1] & dd < bin[2]
        if (sum(s) < 20) next
        res[[length(res) + 1]] <- data.table(family = fam, qtl_chr = ch, group = grp,
          bin = paste0(bin[1], "-", bin[2]), n_pairs = sum(s), mean_r2 = mean(rr[s]))
      }
    }
    ## block fragmentation (genome-wide) for this group
    for (li in idx) {
      g <- G[, li]
      chr_vec <- mk$chr; pos_vec <- mk$pos
      ok <- !is.na(g)
      if (sum(ok) < 50) next
      runs <- rle(g[ok]); chrs <- chr_vec[ok]
      ## count ancestry switches only within chromosome
      sw <- sum(runs$lengths > 1)  # approximate number of runs - 1 per chromosome
      blk[[length(blk) + 1]] <- data.table(family = fam, qtl_chr = ch, group = grp,
        line = colnames(G)[li], n_runs = length(runs$lengths), mean_run_markers = mean(runs$lengths))
    }
  }
}
R <- rbindlist(res); B <- rbindlist(blk)
fwrite(R, file.path(TAB, "Table67_LD_by_QTL_genotype.csv"))
fwrite(B, file.path(TAB, "Table68_blocks_by_QTL_genotype.csv"))
cat("\nLD means by group/bin:\n")
print(R[, .(mean_r2 = round(mean(mean_r2), 3), n = .N), by = .(group, bin)][order(bin, group)])
cat("\nblock stats by group:\n")
print(B[, .(n_runs = round(mean(n_runs), 1), mean_run = round(mean(mean_run_markers), 1), n = .N), by = group])
## significance tests
R[, grp := ifelse(group == "carrier", "carrier", "non")]
W <- dcast(R[, .(mean_r2 = mean(mean_r2)), by = .(family, bin, grp)], family + bin ~ grp, value.var = "mean_r2")
W <- W[!is.na(carrier) & !is.na(non)]
cat("\npaired Wilcoxon across family x bin (carrier vs non): p =",
    signif(wilcox.test(W$carrier, W$non, paired = TRUE)$p.value, 3), " n =", nrow(W), "\n")
cat("mean difference (carrier - non):", round(mean(W$carrier - W$non), 4), "\n")
Bm <- B[, .(n_runs = mean(n_runs), mean_run = mean(mean_run_markers)), by = .(family, group)]
Bw <- dcast(Bm, family ~ group, value.var = "n_runs")
Bw <- Bw[!is.na(carrier) & !is.na(non)]
cat("block number paired Wilcoxon p =", signif(wilcox.test(Bw$carrier, Bw$non, paired = TRUE)$p.value, 3), "\n")
Bw2 <- dcast(Bm, family ~ group, value.var = "mean_run")
Bw2 <- Bw2[!is.na(carrier) & !is.na(non)]
cat("block length paired Wilcoxon p =", signif(wilcox.test(Bw2$carrier, Bw2$non, paired = TRUE)$p.value, 3), "\n")
