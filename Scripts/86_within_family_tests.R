## 86_within_family_tests.R -- strictly within-family tests of LD and block fragmentation
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
norm <- function(x) gsub("[.-][0-9][ABD]$", "", x)
read_geno <- function(ds) {
  fp <- file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")); if (!file.exists(fp)) return(NULL)
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]; G
}
r2pair <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) < 15) return(NA_real_); cor(x[ok], y[ok])^2 }
DIR <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
LD <- list(); BL <- list(); tests <- list()
for (i in seq_len(nrow(DIR))) {
  fam <- DIR$family[i]; ch <- DIR$chr[i]
  ds_all <- unique(MK[family == fam & chr == ch]$dataset)
  for (cand in ds_all) {
    G <- read_geno(cand); if (is.null(G)) next
    dd <- MK[dataset == cand & chr == ch]; dd <- dd[marker %in% rownames(G)]
    if (nrow(dd) < 30) next
    dd <- dd[order(cM)]
    pk <- dd[which.min(abs(cM - DIR$pos[i]))]$marker
    gt <- G[pk, ]
    highA <- DIR$direction[i] != "Watkins"
    car <- if (highA) gt == 0 else gt == 1
    non <- if (highA) gt == 1 else gt == 0
    if (sum(car, na.rm = TRUE) < 10 || sum(non, na.rm = TRUE) < 10) next
    ## LD within family on QTL chromosome
    pos <- dd$ppos / 1e6
    for (g in c("carrier", "non")) {
      idx <- if (g == "carrier") which(car) else which(non)
      sub <- G[dd$marker, idx, drop = FALSE]
      set.seed(1); np <- 500
      a <- sample(seq_along(dd$marker), np, TRUE); b <- sample(seq_along(dd$marker), np, TRUE)
      k <- a != b; a <- a[k]; b <- b[k]
      dist <- abs(pos[a] - pos[b])
      rr <- mapply(function(u, v) r2pair(sub[u, ], sub[v, ]), a, b)
      for (bin in list(c(0,5), c(5,10), c(10,20))) {
        s2 <- !is.na(rr) & dist >= bin[1] & dist < bin[2]
        if (sum(s2) < 20) next
        LD[[length(LD) + 1]] <- data.table(family = fam, group = g, bin = paste0(bin[1], "-", bin[2]),
          mean_r2 = mean(rr[s2]), n = sum(s2))
      }
    }
    ## blocks on QTL chromosome (runs of >=2 markers)
    loc <- list()
    for (li in seq_len(ncol(G))) {
      gg <- G[dd$marker, li]; ok <- !is.na(gg)
      if (sum(ok) < 20) next
      r <- rle(gg[ok])
      loc[[length(loc) + 1]] <- data.table(family = fam, line = colnames(G)[li],
        group = ifelse(li %in% which(car), "carrier", ifelse(li %in% which(non), "non", NA)),
        n_runs = sum(r$lengths >= 2), mean_run = mean(r$lengths[r$lengths >= 2]))
    }
    ## within-family tests
    B <- rbindlist(loc)
    BL[[length(BL) + 1]] <- B
    B <- B[!is.na(group)]
    if (nrow(B) > 10) {
      p1 <- tryCatch(wilcox.test(n_runs ~ group, data = B)$p.value, error = function(e) NA)
      p2 <- tryCatch(wilcox.test(mean_run ~ group, data = B)$p.value, error = function(e) NA)
      tests[[length(tests) + 1]] <- data.table(family = fam, chr = ch, trait = DIR$trait[i],
        n_car = sum(car, na.rm = TRUE), n_non = sum(non, na.rm = TRUE),
        runs_car = mean(B[group == "carrier"]$n_runs), runs_non = mean(B[group == "non"]$n_runs),
        len_car = mean(B[group == "carrier"]$mean_run), len_non = mean(B[group == "non"]$mean_run),
        p_runs = signif(p1, 3), p_len = signif(p2, 3))
    }
    break
  }
}
L2 <- rbindlist(LD); B2 <- rbindlist(BL); T2 <- rbindlist(tests)
fwrite(L2, file.path(TAB, "Table76_within_family_LD.csv"))
fwrite(B2, file.path(TAB, "Table77_within_family_blocks.csv"))
fwrite(T2, file.path(TAB, "Table78_within_family_tests.csv"))
cat("families tested:", nrow(T2), "\n")
print(T2)
cat("\nLD within family (carrier vs non):\n")
print(L2[, .(r2 = round(mean(mean_r2), 3)), by = .(bin, group)][order(bin, group)])