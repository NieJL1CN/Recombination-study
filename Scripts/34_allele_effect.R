## 34_allele_effect.R -- allele classes, segregation ratio and CO effects
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])
EV <- fread(file.path(TAB, "Table18c_CO_events.csv"))
CO <- EV[, .(co = .N), by = .(dataset, id)]
fams <- c("ParW750_WATDE0108_watseq_axiom", "RobW181_WATDE0022_dfw_axiom", "ParW007_WATDE0002_watseq_axiom")
all <- list()
for (f in fams) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  raw <- as.matrix(g[, -1, drop = FALSE])
  G <- t(apply(raw, 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 2; o }))
  Hs <- t(apply(raw, 2, function(z) trimws(z) == "H"))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  rownames(Hs) <- colnames(x)[-1]; colnames(Hs) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  Hs <- Hs[match(mk$marker, rownames(Hs)), , drop = FALSE]
  d <- CO[dataset == f]; d <- d[match(colnames(G), id)]
  D <- data.table(dataset = f, marker = mk$marker, chr = mk$chr, ppos = mk$ppos,
                  n0 = rowSums(G == 0, na.rm = TRUE),
                  n2 = rowSums(G == 2, na.rm = TRUE),
                  nH = rowSums(Hs, na.rm = TRUE),
                  nNA = rowSums(is.na(G)))
  D[, freq2 := n2 / (n0 + n2)]
  D[, chi := pchisq(((n2 - n0)^2) / (n0 + n2), 1, lower.tail = FALSE)]
  eff <- t(apply(G, 1, function(v) {
    c(mean(d$co[v == 0], na.rm = TRUE), mean(d$co[v == 2], na.rm = TRUE)) }))
  D[, mean0 := eff[, 1]]
  D[, mean2 := eff[, 2]]
  D[, diff := mean2 - mean0]
  D[, family := sub("_WATDE.*$", "", f)]
  all[[f]] <- D
  cat(f, " markers", nrow(D), " median freq2", round(median(D$freq2, na.rm = TRUE), 3),
      " distorted(p<0.01)", round(mean(D$chi < 0.01, na.rm = TRUE), 3),
      " mean_H", round(mean(D$nH), 2), "\n")
}
A <- rbindlist(all)
fwrite(A, file.path(TAB, "Table21_allele_effects.csv"))
p1 <- ggplot(A, aes(x = freq2)) + geom_histogram(bins = 40, fill = "#2C6E9B", colour = NA) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ family, nrow = 1) + theme_classic(base_size = 7) +
  labs(x = "allele frequency (BB)", y = "markers") + coord_cartesian(xlim = c(0, 1))
ggsave(file.path(FIG, "Fig4a_allele_frequency.png"), p1, width = 6.4, height = 2.4, dpi = 600)
lim <- max(abs(A$diff), na.rm = TRUE)
p2 <- ggplot(A, aes(x = freq2, y = diff)) + geom_point(size = 0.5, alpha = 0.4, colour = "#2C6E9B") +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
  facet_wrap(~ family, nrow = 1) + theme_classic(base_size = 7) +
  labs(x = "allele frequency (BB)", y = "CO effect (BB - AA)") + ylim(-lim, lim)
ggsave(file.path(FIG, "Fig4b_effect_vs_freq.png"), p2, width = 6.4, height = 2.4, dpi = 600)
print(A[, .(n = .N, med_eff = round(median(abs(diff), na.rm = TRUE), 2)),
          by = .(family, freqbin = cut(freq2, c(0, .3, .4, .6, .7, 1)))][order(family, freqbin)])
