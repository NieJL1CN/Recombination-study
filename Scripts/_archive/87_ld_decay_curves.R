## 87_ld_decay_curves.R -- LD decay curves and threshold distances (carrier vs non)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
pair_r2 <- function(M) {
  n <- nrow(M); p <- ncol(M)
  if (p < 5) return(NULL)
  set.seed(1); np <- min(20000, p * (p - 1) / 2)
  a <- sample(p, np, TRUE); b <- sample(p, np, TRUE); k <- a != b; a <- a[k]; b <- b[k]
  out <- numeric(length(a))
  for (i in seq_along(a)) {
    x <- M[, a[i]]; y <- M[, b[i]]; ok <- !is.na(x) & !is.na(y)
    out[i] <- if (sum(ok) >= 15) cor(x[ok], y[ok])^2 else NA_real_
  }
  data.table(i = a, j = b, r2 = out)
}
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
res <- list(); curves <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    G <- read_geno(cand); if (is.null(G)) next
    dd <- MK[dataset == cand & chr == ch]; dd <- dd[marker %in% rownames(G)]
    if (nrow(dd) < 40) next
    dd <- dd[order(cM)]
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    gt <- G[pk, ]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == 0 else gt == 1
    non <- if (highA) gt == 1 else gt == 0
    if (sum(car, na.rm = TRUE) < 12 || sum(non, na.rm = TRUE) < 12) next
    dist <- as.matrix(dist(dd$ppos / 1e6))
    for (g in c("carrier", "non")) {
      idx <- if (g == "carrier") which(car) else which(non)
      M <- t(G[dd$marker, idx, drop = FALSE])
      pr <- pair_r2(M); if (is.null(pr)) next
      pr[, d := dist[cbind(i, j)]]
      pr <- pr[!is.na(r2) & d > 0.05 & d <= 50]
      if (nrow(pr) < 200) next
      curves[[length(curves) + 1]] <- data.table(family = fam, group = g, d = pr$d, r2 = pr$r2)
      fit <- try(loess(r2 ~ d, data = pr, span = 0.5), silent = TRUE)
      if (inherits(fit, "try-error")) next
      grid <- seq(0.1, 50, by = 0.1)
      pred <- predict(fit, grid)
      thr <- function(th) { k <- which(pred <= th)[1]; ifelse(is.na(k), NA_real_, grid[k]) }
      res[[length(res) + 1]] <- data.table(family = fam, chr = ch, group = g,
        n_lines = length(idx), n_pairs = nrow(pr),
        d50 = thr(0.5), d30 = thr(0.3), d20 = thr(0.2))
    }
    break
  }
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table79_LD_decay_distances.csv"))
print(R)
C <- rbindlist(curves)
C[, grp := factor(group, levels = c("non", "carrier"))]
p <- ggplot(C[, .(r2 = mean(r2)), by = .(group, bin = cut(d, seq(0, 50, 1), labels = FALSE))], aes(bin, r2, colour = group)) +
  geom_line(linewidth = 0.4) + theme_classic(base_size = 7) +
  labs(x = "physical distance (Mb)", y = expression(r^2), colour = NULL) +
  theme(legend.position = c(0.8, 0.8))
ggsave(file.path(FIG, "FigS_LD_decay_carrier_vs_non.png"), p, width = 3.4, height = 2.4, dpi = 600)
cat("done\n")
