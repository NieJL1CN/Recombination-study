## ============================================================================
## 05_qtl_family.R -- per-family QTL with the consensus map vs the original map
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(qtl)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
NPERM <- 500

L <- readRDS(file.path(DER, "maps.rds"))
CM <- L$consensus
PH <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
PH[, family := sub("_.*$", "", line)]

parse_pop <- function(f) {
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  h <- x[1:2, ]
  g <- x[-(1:2), ]
  ids <- g[[1]]
  G <- as.matrix(g[, -1, drop = FALSE])
  G <- t(apply(G, 2, function(z) {
    z <- trimws(z)
    out <- rep(NA_integer_, length(z))
    out[z %in% c("A", "a", "1", "AA")] <- 1L
    out[z %in% c("B", "b", "2", "BB")] <- 2L
    out
  }))
  rownames(G) <- colnames(x)[-1]
  colnames(G) <- ids
  map <- data.table(marker = rownames(G), chr = sub("^([1-7][ABD]).*$", "\\1", h[2, -1]),
                    pos = suppressWarnings(as.numeric(h[3, -1])))
  list(G = G, map = map)
}

scan_one <- function(G, map, y, nperm) {
  map <- map[!is.na(map$chr) & !is.na(map$pos), , drop = FALSE]
  map$chr <- trimws(as.character(map$chr))
  if (nrow(map) < 50) return(list(peak = NULL, thr = NA_real_))
  gn <- list()
  for (ch in unique(map$chr)) {
    sel <- which(map$chr == ch)
    m <- t(G[sel, , drop = FALSE]); colnames(m) <- map$marker[sel]; rownames(m) <- colnames(G)
    gn[[ch]] <- list(data = m, map = setNames(map$pos[sel], map$marker[sel])); class(gn[[ch]]) <- "A"
  }
  cr <- list(geno = gn, pheno = data.frame(y = y, row.names = colnames(G)), covar = NULL,
             chrtype = setNames(rep("A", length(gn)), names(gn)), cross = "riself")
  class(cr) <- c("riself", "cross")
  cr <- calc.genoprob(cr, step = 1, error.prob = 1e-4, map.function = "haldane")
  Dm <- do.call(cbind, lapply(cr$geno, function(z) z$prob[, , 2]))
  Dm[is.na(Dm)] <- 0.5
  chr <- unlist(lapply(names(cr$geno), function(ch) rep(ch, dim(cr$geno[[ch]]$prob)[2])), use.names = FALSE)
  ps <- unlist(lapply(cr$geno, function(z) as.numeric(attr(z$prob, "map"))), use.names = FALSE)
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8 & !is.na(ssx)
  Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]
  chr <- chr[keep]; ps <- ps[keep]
  if (!length(ssx)) return(list(peak = NULL, thr = NA_real_))
  yc <- y - mean(y); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999
  lod <- -n / 2 * log10(1 - r2)
  Y <- matrix(yc[sample.int(n, n * nperm, replace = TRUE)], n, nperm)
  Yc <- sweep(Y, 2, colMeans(Y), "-")
  R2p <- crossprod(Dc, Yc)^2 / (ssx * sum(yc^2)); R2p[R2p > 0.999999] <- 0.999999
  th <- quantile(apply(-n / 2 * log10(1 - R2p), 2, max), 0.95)
  dt <- data.table(chr = chr, pos = ps, lod = lod, r2 = r2, eff = num / ssx)
  if (max(lod) <= th) return(list(peak = NULL, thr = as.numeric(th)))
  pk <- dt[lod > th][, .SD[which.max(lod)], by = chr]
  pk[, threshold := as.numeric(th)]
  list(peak = pk, thr = as.numeric(th))
}

fs <- list.files(GEN, pattern = "_rqtl_genotypes\\.csv$", full.names = TRUE)
res <- list()
for (f in fs) {
  tg <- sub("_rqtl_genotypes\\.csv$", "", basename(f))
  fam <- sub("_.*$", "", tg)
  ph <- PH[family == fam]
  if (nrow(ph) < 60) next
  pp <- parse_pop(f)
  keep <- colnames(pp$G) %in% ph$line
  if (sum(keep) < 60) next
  if (sum(keep) < 60) next
  G0 <- pp$G[, keep, drop = FALSE]
  mp0 <- pp$map
  C <- CM[marker %in% rownames(G0)]
  if (uniqueN(C$panel) == 0) next
  pnl <- C[, .N, by = panel][order(-N)]$panel[1]
  C <- C[panel == pnl, .(marker, chr, pos = cM)]
  C <- C[!duplicated(marker)]
  traits <- c("PH_M_cm", "GY_Calc_tha", "GW_M_g1000grn")
  for (tr in traits) {
    y <- ph[match(colnames(G0), line)][[tr]]
    ok <- !is.na(y)
    if (sum(ok) < 60 || sd(y[ok], na.rm = TRUE) == 0) next
    for (map_type in c("consensus", "original")) {
      if (map_type == "consensus") {
        mm <- C[marker %in% rownames(G0)]
        mm <- mm[order(match(chr, CHR), pos)]
        mm[, pos := pos + (seq_len(.N)) * 1e-4, by = .(chr, pos)]
      } else {
        mm <- mp0[marker %in% rownames(G0)]
        mm[, chr := trimws(as.character(chr))]
        mm <- mm[!is.na(chr) & chr %in% CHR]
        mm <- mm[order(match(chr, CHR), pos)]
      }
      idx <- match(mm$marker, rownames(G0)); okm <- !is.na(idx)
      mm <- mm[okm]; idx <- idx[okm]
      if (length(idx) < 50) next
      G <- G0[idx, , drop = FALSE]
      yk <- y[ok]
      # keep lines with phenotype; drop all-NA marker rows handled by calc.genoprob
      Gk <- G[, ok, drop = FALSE]
      r <- try(scan_one(Gk, as.data.frame(mm), yk, NPERM), silent = TRUE)
      if (inherits(r, "try-error")) { cat("ERR", map_type, tg, tr, as.character(r), "\n"); next }
      if (is.null(r$peak)) next
      r$peak[, dataset := tg][, trait := tr][, map := map_type][, panel := pnl]
      res[[paste(tg, tr, map_type)]] <- r$peak
    }
  }
  cat("done", tg, "\n")
}
Q <- rbindlist(res)
fwrite(Q, file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))
cat("QTL rows:", nrow(Q), "\n")
print(Q[, .(n = .N, med_lod = round(median(lod), 2)), by = .(map, trait)][order(map, trait)])