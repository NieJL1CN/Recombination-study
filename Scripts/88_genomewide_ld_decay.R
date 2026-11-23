## 88_genomewide_ld_decay.R -- genome-wide LD decay for carriers vs non-carriers
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
CM <- L$consensus[, .(marker, chr, panel, pos = cM)][!duplicated(marker)]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
res <- list(); curves <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    G <- read_geno(cand); if (is.null(G)) next
    dd <- MK[dataset == cand]
    dd <- dd[marker %in% rownames(G)]
    if (nrow(dd) < 200) next
    pk <- dd[chr == ch][which.min(abs(cM - DIR$pos[i]))]$marker
    if (!length(pk)) next
    gt <- G[pk, ]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == 0 else gt == 1
    non <- if (highA) gt == 1 else gt == 0
    if (sum(car, na.rm = TRUE) < 12 || sum(non, na.rm = TRUE) < 12) next
    for (g in c("carrier", "non")) {
      idx <- if (g == "carrier") which(car) else which(non)
      Mk <- G[dd$marker, idx, drop = FALSE]
      set.seed(1)
      pr <- list()
      for (cc in CHR) {
        sel <- which(dd$chr == cc); if (length(sel) < 10) next
        p <- length(sel); np <- 2500
        a <- sample(p, np, TRUE); b <- sample(p, np, TRUE); k <- a != b; a <- a[k]; b <- b[k]
        pr[[cc]] <- data.table(a = sel[a], b = sel[b], d = abs(dd$ppos[sel[a]] - dd$ppos[sel[b]]) / 1e6)
      }
      P <- rbindlist(pr)
      P <- P[d > 0.05 & d <= 50]
      P[, r2 := mapply(function(u, v) {
        x <- Mk[u, ]; y <- Mk[v, ]; ok <- !is.na(x) & !is.na(y)
        if (sum(ok) >= 15) cor(x[ok], y[ok])^2 else NA_real_ }, a, b)]
      P <- P[!is.na(r2)]
      if (nrow(P) < 300) next
      curves[[length(curves) + 1]] <- data.table(family = fam, group = g, d = P$d, r2 = P$r2)
      fit <- loess(r2 ~ d, data = P, span = 0.5)
      grid <- seq(0.1, 50, 0.1); pred <- predict(fit, grid)
      thr <- function(th) { k <- which(pred <= th)[1]; ifelse(is.na(k), NA_real_, grid[k]) }
      res[[length(res) + 1]] <- data.table(family = fam, chr = ch, group = g,
        n_lines = length(idx), n_pairs = nrow(P), d50 = thr(0.5), d30 = thr(0.3), d20 = thr(0.2))
    }
    break
  }
}
R <- rbindlist(res); C <- rbindlist(curves)
fwrite(R, file.path(TAB, "Table80_genomewide_LD_decay.csv"))
Cb <- C[, .(r2 = mean(r2), n = .N), by = .(group, bin = cut(d, seq(0, 50, 1), labels = FALSE))]
fwrite(Cb, file.path(TAB, "Table80b_LD_curves_binned.csv"))
print(R)
W <- dcast(R, family ~ group, value.var = "d50"); W <- W[!is.na(carrier) & !is.na(non)]
if (nrow(W) > 3) cat("paired d50 (carrier vs non): p =", signif(wilcox.test(W$carrier, W$non, paired = TRUE)$p.value, 3),
                     " mean diff =", round(mean(W$carrier - W$non), 2), "Mb\n")
C[, bin := cut(d, seq(0, 50, 1), labels = FALSE)]
S <- C[, .(r2 = mean(r2)), by = .(group, bin)]
S[, d := bin]
p <- ggplot(S, aes(d, r2, colour = group)) + geom_line(linewidth = 0.5) +
  theme_classic(base_size = 7) + labs(x = "physical distance (Mb)", y = expression(r^2), colour = NULL) +
  theme(legend.position = c(0.8, 0.8))
ggsave(file.path(FIG, "FigS_Genomewide_LD_carrier_vs_non.png"), p, width = 3.4, height = 2.4, dpi = 600)
cat("done\n")