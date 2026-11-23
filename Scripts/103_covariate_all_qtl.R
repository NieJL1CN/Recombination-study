## 103_covariate_all_qtl.R -- HD covariate test for all 9 significant CO-QTL
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]
CO[, middle := R2a + R2b]
CO[, centromere := C]
CO[, whole := R1 + R2a + C + R2b + R3]
HD <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
P20 <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P20[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
lod_profile <- function(Dm, y) {
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8; Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]
  yc <- y - mean(y); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999
  list(lod = -n/2 * log10(1 - r2), keep = keep)
}
res <- list(); prof7 <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]; tr <- DIR$trait[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  G <- NULL; ds <- NA; best <- 0
  for (cand in ds_all) {
    gg <- read_geno(cand); if (is.null(gg)) next
    n <- sum(MK[dataset == cand & chr == ch]$marker %in% rownames(gg))
    if (n > best) { best <- n; G <- gg; ds <- cand }
  }
  if (is.null(G) || best < 20) next
  dd <- MK[dataset == ds & chr == ch]; dd <- dd[marker %in% rownames(G)][order(cM)]
  if (nrow(dd) < 20) next
  y <- CO[match(colnames(G), id)][[tr]]
  cv <- HD$Hd_dto_day[match(colnames(G), HD$line)]
  if (sum(!is.na(y)) < 40 || sum(!is.na(cv)) < 40) next
  cv[is.na(cv)] <- mean(cv, na.rm = TRUE)
  Dm <- t(G[dd$marker, , drop = FALSE]); Dm[is.na(Dm)] <- 0.5
  ok <- !is.na(y); Dm <- Dm[ok, , drop = FALSE]; y2 <- y[ok]; cv2 <- cv[ok]
  p0 <- lod_profile(Dm, y2)
  XX <- cbind(1, cv2); PX <- XX %*% solve(crossprod(XX), t(XX))
  DmC <- Dm - PX %*% Dm; ys <- y2 - PX %*% y2
  p1 <- lod_profile(DmC, ys)
  pos <- dd$cM[p0$keep]
  k0 <- which.max(p0$lod); k1 <- which.max(p1$lod)
  thr <- P20[family == fam & trait == tr & chr == ch]$thr
  res[[length(res)+1]] <- data.table(family = fam, chr = ch, trait = tr,
    lod_nocov = round(max(p0$lod), 2), lod_cov = round(max(p1$lod), 2),
    delta = round(max(p0$lod) - max(p1$lod), 2),
    thr = if (length(thr)) round(thr[1], 2) else NA_real_,
    pos_nocov = round(pos[k0], 1), pos_cov = round(pos[k1], 1))
  if (fam == "ParW352" && ch == "7D") {
    prof7[[1]] <- data.table(pos = pos, lod = p0$lod, model = "no covariate")
    prof7[[2]] <- data.table(pos = pos, lod = p1$lod, model = "Hd covariate")
  }
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table99_HD_covariate_all_QTL.csv"))
print(R)
P <- rbindlist(prof7)
thr7 <- R[family == "ParW352" & chr == "7D"]$thr
p <- ggplot(P, aes(pos, lod, colour = model)) + geom_line(linewidth = 0.4) +
  geom_hline(yintercept = thr7, linetype = 2, linewidth = 0.25, colour = "grey30") +
  theme_classic(base_size = 7) + labs(x = "7D position (cM)", y = "LOD", colour = NULL) +
  theme(legend.position = c(0.75, 0.85))
ggsave(file.path(FIG, "FigS_7D_LOD_Hd_covariate.png"), p, width = 3.6, height = 2.4, dpi = 600)
cat("done\n")