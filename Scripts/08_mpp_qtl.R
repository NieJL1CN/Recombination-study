## ============================================================================
## 08_mpp_qtl.R -- multi-parent QTL with statgenMPP
##   populations = families of one marker panel, parents coded A/B
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(statgenMPP)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
WORK <- file.path(ROOT, "Source/Derived/mpp"); dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
TRAIT <- "GY_Calc_tha"
PANEL <- "AXIOM.1"

L <- readRDS(file.path(DER, "maps.rds")); CM <- L$consensus
PH <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv")); PH[, family := sub("_.*$", "", line)]
CP <- CM[panel == PANEL, .(marker, chr, pos = cM)][!duplicated(marker)]

fs <- list.files(GEN, pattern = "_rqtl_genotypes\\.csv$", full.names = TRUE)
pops <- character(0); genos <- list(); phenos <- list()
for (f in fs) {
  tg <- sub("_rqtl_genotypes\\.csv$", "", basename(f)); fam <- sub("_.*$", "", tg)
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  h <- x[1:2, ]; g <- x[-(1:2), ]
  ids <- g[[1]]; mk <- colnames(x)[-1]
  keep <- mk %in% CP$marker
  if (sum(keep) < 100) next
  M <- matrix(NA_integer_, nrow(g), sum(keep))
  sub <- as.matrix(g[, -1, drop = FALSE][, keep, drop = FALSE])
  for (j in seq_len(ncol(sub))) {
    z <- trimws(sub[, j]); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L
    M[, j] <- o
  }
  colnames(M) <- mk[keep]; rownames(M) <- ids
  ph <- PH[family == fam]
  if (nrow(ph) < 60) next
  pr <- ph[match(ids, line)]
  ok <- !is.na(pr[[TRAIT]]) & rowSums(!is.na(M)) > 50
  if (sum(ok) < 60) next
  pops <- c(pops, tg)
  genos[[tg]] <- M[ok, , drop = FALSE]
  phenos[[tg]] <- data.table(genotype = ids[ok], trait = pr[[TRAIT]][ok])
  cat("pop", tg, "lines", sum(ok), "markers", ncol(M), "\n")
}
if (length(pops) < 2) stop("need >=2 populations")

common <- Reduce(intersect, lapply(genos, colnames))
cat("common markers:", length(common), "
")
genos <- lapply(genos, function(M) M[, common, drop = FALSE])
allm <- sort(common)
mp <- CP[marker %in% allm]
mp <- mp[!is.na(chr) & !is.na(pos)]
mp <- unique(mp, by = "marker")
mp <- mp[order(match(chr, CHR), pos)]
mp <- mp[, if (.N >= 2) .SD, by = chr]
mp[, pos := pos + seq_len(.N) * 1e-4, by = chr]
setcolorder(mp, c("marker", "chr", "pos"))
fwrite(mp, file.path(WORK, "mapfile.txt"), sep = "\t", col.names = FALSE)

files <- character(0)
for (p in pops) {
  M <- genos[[p]]
  D <- data.table(ID = c("A", "B", rownames(M)))
  G <- rbind(rep(1L, ncol(M)), rep(2L, ncol(M)), M)
  colnames(G) <- colnames(M)
  D <- cbind(D, as.data.table(G))
  fp <- file.path(WORK, paste0(p, ".txt"))
  fwrite(D, fp, sep = "\t", quote = FALSE, na = "NA")
  files <- c(files, fp)
}
PHN <- rbindlist(phenos)
PHN[, genotype := paste0(genotype)]
setnames(PHN, "trait", TRAIT)
fwrite(PHN, file.path(WORK, "pheno.txt"), sep = "\t")

res <- try(calcIBDMPP(crossNames = pops, markerFiles = files, pheno = as.data.frame(PHN),
                      popType = "F6", mapFile = file.path(WORK, "mapfile.txt"),
                      evalDist = 10, grid = TRUE, verbose = TRUE), silent = TRUE)
if (inherits(res, "try-error")) {
  cat("calcIBDMPP ERROR:", as.character(res), "\n")
} else {
  saveRDS(res, file.path(DER, "mpp_ibd.rds"))
  cat("IBD computed\n"); print(summary(res))
  q <- try(selQTLMPP(res, trait = TRAIT, QTLwindow = 10, threshold = 3, verbose = FALSE), silent = TRUE)
  if (inherits(q, "try-error")) cat("selQTLMPP ERROR:", as.character(q), "\n") else {
    saveRDS(q, file.path(DER, "mpp_qtl.rds")); print(q)
  }
}