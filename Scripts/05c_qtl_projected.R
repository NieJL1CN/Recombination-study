## ============================================================================
## 05c_qtl_projected.R -- identical marker set in both maps:
##   family markers projected onto consensus coordinates (isotonic mapping)
##   vs the same markers on the family's own map
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(qtl)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
NPERM <- 500
L <- readRDS(file.path(DER, "maps_min1.rds")); CM <- L$consensus
PH <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv")); PH[, family := sub("_.*$", "", line)]

scan2 <- function(G, map, y, nperm = NPERM) {
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
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8; Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]; chr <- chr[keep]; ps <- ps[keep]
  yc <- y - mean(y); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999; lod <- -n / 2 * log10(1 - r2)
  Y <- matrix(yc[sample.int(n, n * nperm, TRUE)], n, nperm)
  Yc <- sweep(Y, 2, colMeans(Y), "-")
  R2p <- crossprod(Dc, Yc)^2 / (ssx * sum(yc^2)); R2p[R2p > 0.999999] <- 0.999999
  th <- quantile(apply(-n / 2 * log10(1 - R2p), 2, max), 0.95)
  data.table(chr = chr, pos = ps, lod = lod, r2 = r2, eff = num / ssx, thr = as.numeric(th))
}

fs <- list.files(GEN, pattern = "_rqtl_genotypes\\.csv$", full.names = TRUE)
out <- list()
for (f in fs) {
  tg <- sub("_rqtl_genotypes\\.csv$", "", basename(f)); fam <- sub("_.*$", "", tg)
  ph <- PH[family == fam]; if (nrow(ph) < 60) next
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  h <- x[1:2, ]; g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mp <- data.table(marker = rownames(G),
                   chr = trimws(sub("^([1-7][ABD]).*$", "\\1", unlist(h[1, -1]))),
                   own = suppressWarnings(as.numeric(unlist(h[2, -1]))))
  mp <- mp[chr %in% CHR & !is.na(own)]
  C <- CM[marker %in% mp$marker]; if (!nrow(C)) next
  pnl <- C[, .N, by = panel][order(-N)]$panel[1]
  C <- C[panel == pnl, .(marker, cns = cM)]
  mp <- merge(mp, C, by = "marker", all.x = TRUE)
  sh <- mp[!is.na(cns)][order(own)]
  if (nrow(sh) < 30) next
  iso <- isoreg(sh$own, sh$cns)
  f_map <- approxfun(x = iso$x, y = iso$y, rule = 2, ties = "ordered")
  mp[, proj := f_map(own)]
  mp[is.na(proj), proj := own]
  keep <- colnames(G) %in% ph$line; if (sum(keep) < 60) next
  G0 <- G[, keep, drop = FALSE]
  for (tr in c("PH_M_cm", "GY_Calc_tha", "GW_M_g1000grn")) {
    y <- ph[match(colnames(G0), line)][[tr]]; ok <- !is.na(y)
    if (sum(ok) < 60 || sd(y[ok]) == 0) next
    yk <- y[ok]
    for (mt in c("consensus", "original")) {
      mm <- copy(mp)
      mm[, pos := if (mt == "consensus") proj else own]
      mm <- mm[order(match(chr, CHR), pos)]
      mm[, pos := pos + seq_len(.N) * 1e-5, by = chr]
      r <- try(scan2(G0[, ok, drop = FALSE], mm, yk), silent = TRUE)
      if (inherits(r, "try-error")) { cat("ERR", tg, tr, mt, "\n"); next }
      r[, dataset := tg]; r[, trait := tr]; r[, map := mt]; r[, panel := pnl]; r[, n_marker := nrow(mm)]
      out[[paste(tg, tr, mt)]] <- r
    }
  }
  cat("done", tg, "\n")
}
PR <- rbindlist(out)
fwrite(PR, file.path(TAB, "Table2d_QTL_profiles.csv"))
PK <- PR[lod > thr, .SD[which.max(lod)], by = .(dataset, trait, map, chr)]
fwrite(PK, file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))
cat("profiles:", nrow(PR), " peaks:", nrow(PK), "\n")
print(PK[, .(n = .N, med_lod = round(median(lod), 2)), by = .(map, trait)][order(map, trait)])
