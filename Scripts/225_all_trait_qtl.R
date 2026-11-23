## 225_all_trait_qtl.R -- single-marker QTL scan for all available traits
suppressPackageStartupMessages({library(data.table); library(qtl)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
CM <- readRDS(file.path(DER, "maps_min1.rds"))$consensus[, .(marker, chr, panel, pos = cM)]
CM <- CM[!duplicated(marker)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
p43 <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
p2a <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
PH <- merge(p43[, .(line, Hd_dto_day, PH_M_cm, GY_Calc_tha, GW_M_g1000grn)],
            p2a[, .(line, GN_Calc_Grnm2, GWid_M_mm, GLng_M_mm)], by = "line", all = TRUE)
traits <- c("Hd_dto_day","PH_M_cm","GY_Calc_tha","GW_M_g1000grn","GN_Calc_Grnm2","GWid_M_mm","GLng_M_mm")
fams <- unique(sub("_WATDE.*$", "", CO$dataset))
NPERM <- 200
out <- list()
for (fam in fams) {
  ds <- unique(CO[grepl(paste0("^", fam, "_"), dataset)]$dataset)
  for (d in ds) {
    fp <- file.path(GEN, paste0(d, "_rqtl_genotypes.csv"))
    if (!file.exists(fp)) next
    x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
    Gm <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
      z <- trimws(z); o <- rep(NA_real_, length(z))
      o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
    rownames(Gm) <- colnames(x)[-1]; colnames(Gm) <- g[[1]]
    C <- CM[marker %in% rownames(Gm)]
    if (!nrow(C)) next
    pnl <- C[, .N, by = panel][order(-N)]$panel[1]
    C <- C[panel == pnl]
    Gc <- Gm[C$marker, , drop = FALSE]
    ids <- colnames(Gc)
    Dm <- t(Gc); Dm[is.na(Dm)] <- 0.5
    Dc <- sweep(Dm, 2, colMeans(Dm), "-")
    ssx <- colSums(Dc^2); keep <- ssx > 1e-8
    Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]
    chr <- C$chr[keep]; ps <- C$pos[keep]
    for (tr in traits) {
      y <- PH[[tr]][match(ids, PH$line)]
      if (sum(!is.na(y)) < 40) next
      ok <- !is.na(y); yc <- y[ok] - mean(y[ok]); Do <- Dc[ok, , drop = FALSE]
      n <- length(yc); num <- as.numeric(crossprod(Do, yc))
      r2 <- num^2 / (ssx * sum(yc^2)); r2[r2 > 0.999] <- 0.999
      lod <- -n/2 * log10(1 - r2)
      Y <- matrix(yc[sample.int(n, n*NPERM, TRUE)], n, NPERM); Yc <- sweep(Y, 2, colMeans(Y), "-")
      Rp <- crossprod(Do, Yc)^2 / (ssx * sum(yc^2)); Rp[Rp > 0.999] <- 0.999
      thr <- quantile(apply(-n/2*log10(1-Rp), 2, max), 0.95)
      dt <- data.table(family = fam, dataset = d, trait = tr, chr = chr, pos = ps, lod = lod, r2 = r2, thr = as.numeric(thr))
      dt <- dt[lod >= thr]
      if (nrow(dt)) {
        dt <- dt[order(-lod)][, head(.SD, 1), by = chr]
        out[[length(out)+1]] <- dt
      }
    }
  }
  cat(fam, "done\n")
}
R <- rbindlist(out, fill = TRUE)
fwrite(R, file.path(TAB, "Table175_all_trait_QTL.csv"))
cat("QTL rows", nrow(R), "\n"); print(R[, .N, by = trait])
