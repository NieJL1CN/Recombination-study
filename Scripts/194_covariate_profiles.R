## 194_covariate_profiles.R -- LOD profiles with/without HD for the three families
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(qtl); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
NPERM <- 5
CM <- readRDS(file.path(DER, "maps_min1.rds"))$consensus[, .(marker, chr, panel, pos = cM)]
CM <- CM[!duplicated(marker)]
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]; CO[, middle := R2a + R2b]; CO[, centromere := C]; CO[, whole := R1 + R2a + C + R2b + R3]
HD <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
P20 <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv")); P20[, family := sub("_WATDE.*$", "", dataset)]
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
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
  cr <- qtl::calc.genoprob(cr, step = 1, error.prob = 1e-4, map.function = "haldane")
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
profs <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]; tr <- DIR$trait[i]
  ds_all <- unique(P20[family == fam & chr == ch]$dataset)
  best <- 0; G <- NULL
  for (cand in ds_all) {
    fp <- file.path(GEN, paste0(cand, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
    x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
    Gm <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
      z <- trimws(z); o <- rep(NA_real_, length(z))
      o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
    rownames(Gm) <- colnames(x)[-1]; colnames(Gm) <- g[[1]]
    C <- CM[marker %in% rownames(Gm)]
    if (!nrow(C)) next
    pnl <- C[, .N, by = panel][order(-N)]$panel[1]
    C2 <- C[panel == pnl & chr == ch]
    if (nrow(C2) > best) { best <- nrow(C2); G <- Gm }
  }
  C <- CM[marker %in% rownames(G) & chr == ch]
  Gc <- G[C$marker, , drop = FALSE]
  ids <- colnames(Gc)
  co <- CO[match(ids, id)]
  y <- if (tr == "distal") co$distal else if (tr == "whole") co$whole else if (tr == "centromere") co$centromere else co$middle
  hd <- HD$Hd_dto_day[match(ids, HD$line)]
  p0 <- scan2(Gc, C, y)
  p1 <- scan2(Gc, C, y, cv = hd)
  profs[[i]] <- rbind(data.table(family = fam, chr = ch, pos = p0$pos, lod = p0$lod, model = "no covariate"),
                      data.table(family = fam, chr = ch, pos = p1$pos, lod = p1$lod, model = "Hd covariate"))
  cat(fam, ch, "max:", round(max(p0$lod),2), round(max(p1$lod),2), "\n")
}
P <- rbindlist(profs)
fwrite(P, file.path(TAB, "Table168_HD_profiles.csv"))
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  plot.tag = element_text(size = 9, face = "bold"), legend.text = element_text(size = 6.5), legend.title = element_blank()))
p <- ggplot(P, aes(pos, lod, colour = model)) + geom_line(linewidth = 0.5) +
  facet_wrap(~ family + chr, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("no covariate" = "#1B9E9E", "Hd covariate" = "#E0726A")) +
  labs(x = "Chromosome position (cM)", y = "LOD") + theme(legend.position = "top")
ggsave(file.path(OUT, "panels", "Fig3b_profiles.png"), p, width = 150, height = 60, units = "mm", dpi = 600)
cat("saved\n")
