## 32_co_qtl.R -- QTL for four CO phenotypes (distal, middle, centromere, whole)
suppressPackageStartupMessages({library(data.table); library(qtl); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
NPERM <- 500
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
CM <- readRDS(file.path(DER, "maps_min1.rds"))$consensus[, .(marker, chr, panel, pos = cM)]
CM <- CM[!duplicated(marker)]
RG <- fread(file.path(TAB, paste0("Table22c_CO_rate_individual", Sys.getenv("SUFFIX", ""), ".csv")))
NW <- fread(file.path(TAB, "Table23_informative_windows.csv"))
PH <- RG[, .(dataset, id,
             distal = R1 + R3,
             middle = R2a + R2b,
             centromere = C,
             whole = R1 + R2a + C + R2b + R3)]
PH <- merge(PH, NW, by = c("dataset", "id"), all.x = TRUE)

scan2 <- function(G, map, y, cv = NULL, nperm = NPERM) {
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
  if (!is.null(cv) && sum(!is.na(cv)) > 10) {
    XX <- cbind(1, cv); PX <- XX %*% solve(crossprod(XX), t(XX))
    Dm <- Dm - PX %*% Dm
    ys <- y - PX %*% y
  } else ys <- y
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8; Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]; chr <- chr[keep]; ps <- ps[keep]
  yc <- as.numeric(ys) - mean(ys); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999; lod <- -n / 2 * log10(1 - r2)
  Y <- matrix(yc[sample.int(n, n * nperm, TRUE)], n, nperm); Yc <- sweep(Y, 2, colMeans(Y), "-")
  R2p <- crossprod(Dc, Yc)^2 / (ssx * sum(yc^2)); R2p[R2p > 0.999999] <- 0.999999
  th <- quantile(apply(-n / 2 * log10(1 - R2p), 2, max), 0.95)
  data.table(chr = chr, pos = ps, lod = lod, r2 = r2, eff = num / ssx, thr = as.numeric(th))
}

res <- list(); prof <- list()
for (f in COH) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  C <- CM[marker %in% rownames(G)]
  if (!nrow(C)) next
  pnl <- C[, .N, by = panel][order(-N)]$panel[1]
  C <- C[panel == pnl]
  mk <- C[order(match(chr, CHR), pos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  ph <- PH[dataset == f]; ph <- ph[match(colnames(G), id)]
  for (tr in c("distal", "middle", "centromere", "whole")) {
    y <- ph[[tr]]
    if (sum(!is.na(y)) < 60 || sd(y, na.rm = TRUE) == 0) next
    r <- try(scan2(G, mk, y, cv = ph$n_win), silent = TRUE)
    if (inherits(r, "try-error")) next
    r[, dataset := f][, trait := tr]
    prof[[paste(f, tr)]] <- r
    pk <- r[lod > thr][, .SD[which.max(lod)], by = chr]
    if (nrow(pk)) res[[paste(f, tr)]] <- pk
  }
  cat("done", f, "\n")
}
P <- rbindlist(prof); Q <- rbindlist(res)
fwrite(P, file.path(TAB, paste0("Table20_CO_QTL_profiles", Sys.getenv("SUFFIX", ""), ".csv")))
fwrite(Q, file.path(TAB, paste0("Table20b_CO_QTL_peaks", Sys.getenv("SUFFIX", ""), ".csv")))
cat("peaks:", nrow(Q), "\n"); print(Q[, .N, by = trait])