## 45_ibspy_enrichment.R -- IBSpy haplotype diversity at CO-QTL intervals
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
W <- fread(file.path(TAB, "Table33_IBSpy_window_diversity.csv"))
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
META <- fread(file.path(TAB, "Table30_CO_metaQTL.csv"))[n_families >= 2]
P <- merge(PK, META[, .(chr, cluster)], by = c("chr", "cluster"))
P[, phys := NA_real_]
for (i in seq_len(nrow(P))) {
  d <- MK[family == P$family[i] & chr == P$chr[i]]
  if (nrow(d)) P$phys[i] <- d[which.min(abs(cM - P$pos[i])), ppos]
}
PI <- P[!is.na(phys), .(Mb_start = min(phys) / 1e6, Mb_end = max(phys) / 1e6, n_peaks = .N),
        by = .(chr, cluster)]
PI[, Mb_start := pmax(0, Mb_start - 5)]
PI[, Mb_end := Mb_end + 5]
fwrite(PI, file.path(TAB, "Table34_CO_metaQTL_physical.csv"))
res <- list()
for (i in seq_len(nrow(PI))) {
  ch <- PI$chr[i]; s0 <- PI$Mb_start[i] * 1e6; e0 <- PI$Mb_end[i] * 1e6
  w <- W[chr == ch]
  inq <- w[start >= s0 & end <= e0]
  bg <- w[!(start >= s0 & end <= e0)]
  if (nrow(inq) < 2 || nrow(bg) < 10) next
  wt <- wilcox.test(inq$shannon, bg$shannon)
  res[[i]] <- data.table(chr = ch, cluster = PI$cluster[i],
                         Mb_start = round(s0 / 1e6, 1), Mb_end = round(e0 / 1e6, 1),
                         n_win_in = nrow(inq), n_bg = nrow(bg),
                         shannon_in = round(mean(inq$shannon), 3),
                         shannon_bg = round(mean(bg$shannon), 3),
                         groups_in = round(mean(inq$n_group), 1),
                         groups_bg = round(mean(bg$n_group), 1),
                         p = wt$p.value)
}
R <- rbindlist(res)
R[, p_adj := p.adjust(p, "BH")]
fwrite(R, file.path(TAB, "Table34b_IBSpy_diversity_test.csv"))
print(R)
W[, chr := factor(chr, levels = CHR)]
PI[, chr := factor(chr, levels = CHR)]
p <- ggplot(W, aes(x = start / 1e6, y = shannon)) +
  geom_line(colour = "grey65", linewidth = 0.12) +
  geom_rect(data = PI, aes(xmin = Mb_start, xmax = Mb_end, ymin = -Inf, ymax = Inf),
            fill = "#E08A2E", alpha = 0.2, inherit.aes = FALSE) +
  facet_wrap(~ chr, ncol = 7, scales = "free_x") +
  theme_classic(base_size = 6) +
  theme(strip.text = element_text(size = 5), axis.text = element_text(size = 3.5)) +
  labs(x = "physical position (Mb)", y = "haplotype diversity (Shannon)")
ggsave(file.path(FIG, "Fig4d_IBSpy_diversity_CO_QTL.png"), p, width = 7.4, height = 3.4, dpi = 600)
cat("done\n")
