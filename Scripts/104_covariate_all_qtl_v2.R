## 104_covariate_all_qtl_v2.R -- HD covariate test, same pipeline as CO-QTL scan
suppressPackageStartupMessages({library(data.table); library(qtl); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank)]
MK[, family := sub("_WATDE.*$", "", dataset)]
CM <- L$consensus[, .(marker, chr, panel, pos = cM)][!duplicated(marker)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]
CO[, middle := R2a + R2b]
CO[, centromere := C]
CO[, whole := R1 + R2a + C + R2b + R3]
HD <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
P20 <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P20[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
scan_lod <- function(G, map, y, cv = NULL) {
  gn <- list()
  for (ch in CHR) {
    sel <- which(map$chr == ch); if (!length(sel)) next
    m <- t(G[sel, , drop = FALSE]); colnames(m) <- map$marker[sel]; rownames(m) <- colnames(G)
    gn[[ch]] <- list(data = m, map = setNames(map$pos[sel], map$marker[sel])); class(gn[[ch]]) <- "A"
  }
  cr <- list(geno = gn, pheno = data.frame(y = y, row.names = colnames(G)), covar = NULL,
             chrtype = setNames(rep("A", length(gn)), names(gn)), cross = "riself")
  class(cr) <- c("riself", "cross")
  cr <- calc.genoprob(cr, step = 1, error.prob = 1e-4, map.function = "haldane")
  Dm <- do.call(cbind, lapply(cr$geno, function(z) z$prob[, , 2])); Dm[is.na(Dm)] <- 0.5
  chr <- unlist(lapply(names(cr$geno), function(ch) rep(ch, dim(cr$geno[[ch]]$prob)[2])))
  ps <- unlist(lapply(cr$geno, function(z) as.numeric(attr(z$prob, "map"))))
  if (!is.null(cv)) { XX <- cbind(1, cv); PX <- XX %*% solve(crossprod(XX), t(XX)); Dm <- Dm - PX %*% Dm; ys <- y - PX %*% y } else ys <- y
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8; Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]; chr <- chr[keep]; ps <- ps[keep]
  yc <- as.numeric(ys) - mean(ys); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999
  data.table(chr = chr, pos = ps, lod = -n/2 * log10(1 - r2))
}
res <- list(); profs <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]; tr <- DIR$trait[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  G <- NULL; ds <- NA; best <- 0
  for (cand in ds_all) {
    fp <- file.path(GEN, paste0(cand, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
    x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
    Gm <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
      z <- trimws(z); o <- rep(NA_real_, length(z))
      o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
    rownames(Gm) <- colnames(x)[-1]; colnames(Gm) <- g[[1]]
    n <- sum(MK[dataset == cand & chr == ch]$marker %in% rownames(Gm))
    if (n > best) { best <- n; G <- Gm; ds <- cand }
  }
  if (is.null(G) || best < 20) { cat(fam, ch, "skipped\n"); next }
  C <- CM[marker %in% rownames(G)]
  pnl <- C[, .N, by = panel][order(-N)]$panel[1]; C <- C[panel == pnl]
  mk <- C[order(match(chr, CHR), pos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  y <- CO[match(colnames(G), id)][[tr]]
  cv <- HD$Hd_dto_day[match(colnames(G), HD$line)]
  keep <- !is.na(y) & !is.na(cv)
  if (sum(keep) < 40) { cat(fam, ch, "few samples\n"); next }
  G <- G[, keep, drop = FALSE]; y <- y[keep]; cv <- cv[keep]
  r0 <- scan_lod(G, mk, y); r1 <- scan_lod(G, mk, y, cv)
  d0 <- r0[chr == ch]; d1 <- r1[chr == ch]
  thr <- P20[family == fam & trait == tr & chr == ch]$thr
  res[[length(res)+1]] <- data.table(family = fam, chr = ch, trait = tr,
    lod_nocov = round(max(d0$lod), 2), lod_cov = round(max(d1$lod), 2),
    delta = round(max(d0$lod) - max(d1$lod), 2),
    thr = if (length(thr)) round(thr[1], 2) else NA_real_,
    pos_nocov = round(d0$pos[which.max(d0$lod)], 1), pos_cov = round(d1$pos[which.max(d1$lod)], 1))
  if (fam == "ParW352" && ch == "7D") {
    profs[[1]] <- data.table(pos = d0$pos, lod = d0$lod, model = "no covariate")
    profs[[2]] <- data.table(pos = d1$pos, lod = d1$lod, model = "Hd covariate")
  }
  cat(fam, ch, "done\n")
}
R <- rbindlist(res); fwrite(R, file.path(TAB, "Table99_HD_covariate_all_QTL.csv"))
print(R)
if (length(profs)) {
  P <- rbindlist(profs); thr7 <- R[family == "ParW352" & chr == "7D"]$thr
  p <- ggplot(P, aes(pos, lod, colour = model)) + geom_line(linewidth = 0.4) +
    geom_hline(yintercept = thr7, linetype = 2, linewidth = 0.25, colour = "grey30") +
    theme_classic(base_size = 7) + labs(x = "7D position (cM)", y = "LOD", colour = NULL) +
    theme(legend.position = c(0.75, 0.85))
  ggsave(file.path(FIG, "FigS_7D_LOD_Hd_covariate.png"), p, width = 3.6, height = 2.4, dpi = 600)
}
cat("finished\n")
