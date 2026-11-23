## ============================================================================
## 15_co_genome_dist.R -- CO density along the genome (10 Mb bins, region colour)
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- 5.5; STEP <- 2.75
REG <- c("R1" = "#2C6E9B", "R2a" = "#6A4C93", "C" = "#BDBDBD", "R2b" = "#2E8B57", "R3" = "#E08A2E")
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)

M <- readRDS(file.path(DER, "maps_min1.rds"))$markers
M <- M[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
ev <- list()
for (f in intersect(unique(M$dataset), COH)) {
  fp <- file.path(GEN, paste0(f, "_rqtl_genotypes.csv")); if (!file.exists(fp)) next
  x <- read.csv(fp, check.names = FALSE, colClasses = "character"); g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_integer_, length(z))
    o[z %in% c("A","a","1","AA")] <- 1L; o[z %in% c("B","b","2","BB")] <- 2L; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- mk[order(match(chr, CHR), ppos)]
  for (ch in unique(mk$chr)) {
    s <- mk[chr == ch]; if (nrow(s) < 10) next
    starts <- seq(0, max(s$ppos), by = STEP * 1e6)
    W <- matrix(NA_integer_, length(starts), ncol(G))
    for (k in seq_along(starts)) {
      rr <- which(s$ppos >= starts[k] & s$ppos < starts[k] + WIN * 1e6)
      if (!length(rr)) next
      a <- G[rr, , drop = FALSE]
      nA <- colSums(a == 1, na.rm = TRUE); nB <- colSums(a == 2, na.rm = TRUE)
      W[k, ] <- ifelse(nA > 0 & nB == 0, 1L, ifelse(nB > 0 & nA == 0, 2L, NA_integer_))
    }
    mid <- starts + WIN * 1e6 / 2
    for (j in seq_len(ncol(G))) {
      v <- W[, j]; keep <- !is.na(v); if (sum(keep) < 2) next
      kk <- which(keep); sw <- kk[-1][v[kk[-1]] != v[kk[-length(kk)]]]
      if (!length(sw)) next
      ev[[length(ev) + 1]] <- data.table(dataset = f, chr = ch, ppos = mid[sw])
    }
  }
  cat("done", f, "\n")
}
E <- rbindlist(ev)
E[, bin := floor(ppos / 10e6)]
E <- merge(E, B, by = "chr")
E[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
      fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
CNT <- E[, .(n = .N), by = .(chr, bin, region)]
NID <- unique(E[, .(dataset, chr)])[, .(n_id = .N), by = chr]
CNT <- merge(CNT, NID, by = "chr")
CNT[, co_per_ind := n / n_id]
fwrite(CNT, file.path(TAB, "Table11_CO_10Mb_bins.csv"))

CH <- unique(CNT[, .(chr, bin)])[order(match(chr, CHR), bin)]
CH[, cum := (match(chr, CHR) - 1), by = chr]
CH[, x := as.numeric(factor(chr, levels = CHR)) - 0.5 + (bin * 10 - min(bin * 10)) / 900 * 0.9, by = chr]
CH <- merge(CH, CNT[, .(chr, bin, region, co_per_ind)], by = c("chr", "bin"))
CH[, chr := factor(chr, levels = CHR)]
p <- ggplot(CH, aes(x = as.numeric(chr), y = co_per_ind, fill = region)) +
  geom_col(width = 0.85, colour = "grey35", linewidth = 0.05) +
  scale_fill_manual(values = REG, name = NULL) +
  scale_x_continuous(breaks = seq_along(CHR), labels = CHR) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4),
        legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 5)) +
  labs(x = NULL, y = "CO per individual (10 Mb bin)")
ggsave(file.path(FIG, "Fig3d_CO_genome_distribution.png"), p, width = 6.8, height = 3, dpi = 600)
cat("bins:", nrow(CNT), "\n")