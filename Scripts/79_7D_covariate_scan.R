## 79_7D_covariate_scan.R -- 7D interval mapping with heading date as covariate (ParW352)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(qtl)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
NPERM <- 500
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
CM <- readRDS(file.path(DER, "maps_min1.rds"))$consensus[, .(marker, chr, panel, pos = cM)]
CM <- CM[!duplicated(marker)]
RG <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
PH <- RG[, .(dataset, id, distal = R1 + R3)]
HD <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
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
    Dm <- Dm - PX %*% Dm; ys <- y - PX %*% y
  } else ys <- y
  Dc <- sweep(Dm, 2, colMeans(Dm), "-"); ssx <- colSums(Dc^2)
  keep <- ssx > 1e-8; Dc <- Dc[, keep, drop = FALSE]; ssx <- ssx[keep]; chr <- chr[keep]; ps <- ps[keep]
  yc <- as.numeric(ys) - mean(ys); n <- length(yc)
  num <- as.numeric(crossprod(Dc, yc)); r2 <- num^2 / (ssx * sum(yc^2))
  r2[r2 > 0.999999] <- 0.999999; lod <- -n / 2 * log10(1 - r2)
  data.table(chr = chr, pos = ps, lod = lod, r2 = r2, eff = num / ssx)
}
f <- "ParW352_WATDE0047_watseq_axiom"
fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv"))
x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
  z <- trimws(z); o <- rep(NA_real_, length(z))
  o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
C <- CM[marker %in% rownames(G)]
pnl <- C[, .N, by = panel][order(-N)]$panel[1]
C <- C[panel == pnl]
mk <- C[order(match(chr, CHR), pos)]
G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
ph <- PH[dataset == f]; ph <- ph[match(colnames(G), id)]
y <- ph$distal
cv <- HD$Hd_dto_day[match(colnames(G), HD$line)]
cv[is.na(cv)] <- mean(cv, na.rm = TRUE)
cat("lines:", sum(!is.na(y)), " with Hd:", sum(!is.na(cv)), "\n")
r0 <- scan2(G, mk, y, cv = NULL)
r1 <- scan2(G, mk, y, cv = cv)
d0 <- r0[chr == "7D"]; d1 <- r1[chr == "7D"]
cat("7D peak without covariate: LOD", round(max(d0$lod), 2), "at", d0$pos[which.max(d0$lod)], "cM\n")
cat("7D peak with Hd covariate: LOD", round(max(d1$lod), 2), "at", d1$pos[which.max(d1$lod)], "cM\n")
cat("LOD drop at original peak:", round(d0$lod[which.max(d0$lod)] - d1$lod[which.min(abs(d1$pos - d0$pos[which.max(d0$lod)]))], 2), "\n")
R <- rbind(d0[, model := "no covariate"], d1[, model := "Hd covariate"])
fwrite(R, file.path(TAB, "Table66_7D_LOD_Hd_covariate.csv"))
THR <- 2.94  ## permutation threshold for ParW352 7D distal (Table20)
p <- ggplot(R, aes(pos, lod, colour = model)) + geom_line(linewidth = 0.3) +
  geom_hline(yintercept = THR, linetype = 2, linewidth = 0.25, colour = "grey30") +
  theme_classic(base_size = 7) + labs(x = "7D position (cM)", y = "LOD", colour = NULL) +
  theme(legend.position = c(0.8, 0.85))
ggsave(file.path(FIG, "FigS_7D_LOD_Hd_covariate.png"), p, width = 3.6, height = 2.4, dpi = 600)
cat("saved figure\n")