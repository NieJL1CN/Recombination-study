## 11_meta_map.R -- vertical chromosome + QTL CIs (1-LOD drop) + family columns
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
CONS <- L$consensus[, .(marker, chr, cM, panel)]
TRCOL <- c(PH = "#2C6E9B", GY = "#2E8B57", GW = "#E08A2E")
B <- data.table(chr = CHR,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6)
POS <- unique(L$markers[!is.na(g_rank), .(marker, ppos)])
D2 <- merge(CONS, POS, by = "marker")
D2[, midp := (B$b2[match(chr, B$chr)] + B$b3[match(chr, B$chr)]) / 2]
CENT <- D2[, .SD[which.min(abs(ppos - midp[1]))], by = chr][, .(chr, cen = cM)]
LEN <- CONS[, .(len = max(cM)), by = chr]
MKPOS <- CONS[, .(chr, marker, cM)]
## QTL peaks and 1-LOD-drop CIs from the saved profiles (same rule as qtl::lodint)
Q <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))[map == "consensus"]
PR <- fread(file.path(TAB, "Table2d_QTL_profiles.csv"))[map == "consensus"]
CI <- PR[, {
  r <- .SD[order(pos)]
  pk <- which.max(r$lod)
  thr <- r$lod[pk] - 1.0
  lo <- pk; while (lo > 1 && r$lod[lo - 1] >= thr) lo <- lo - 1
  hi <- pk; while (hi < .N && r$lod[hi + 1] >= thr) hi <- hi + 1
  list(start = r$pos[lo], end = r$pos[hi], peak = r$pos[pk])
}, by = .(dataset, trait, chr)]
Q <- merge(Q[, .(dataset, trait, chr, tr_short, lod)], CI, by = c("dataset", "trait", "chr"))
Q[, family := sub("_WATDE.*$", "", dataset)]
Q[, trait := factor(tr_short, levels = c("PH","GY","GW"))]
Q <- merge(Q, LEN, by = "chr")
Q[, start := pmax(0, pmin(start, len))][, end := pmax(0, pmin(end, len))]
Q[, peak := pmax(0, pmin(peak, len))]
W <- 0.42; R <- W / 2; maxlod <- Q[, max(lod)]
Q[, lab_y := peak]
Q[, gap := len * 0.045]
setorder(Q, chr, peak)
for (ch in unique(Q$chr)) {
  idx <- which(Q$chr == ch)
  if (length(idx) < 2) next
  g <- Q$gap[idx[1]]
  y <- Q$peak[idx]
  for (k in seq_along(idx)[-1]) {
    if (y[k] - y[k - 1] < g) y[k] <- y[k - 1] + g
  }
  Q$lab_y[idx] <- y
}
bars <- merge(LEN, CENT, by = "chr")
back <- bars[, .(chr, xmin = -R, xmax = R, ymin = 0, ymax = len)]
p <- ggplot() +
  geom_rect(data = back, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = "black", linewidth = 0.25) +
  geom_segment(data = MKPOS, aes(x = -R, xend = R, y = cM, yend = cM),
               colour = "grey55", linewidth = 0.05, alpha = 0.6) +
  geom_rect(data = Q, aes(xmin = -R, xmax = R, ymin = start, ymax = end, fill = trait),
            alpha = 0.65, colour = NA) +
  geom_point(data = CENT, aes(x = 0, y = cen), size = 0.5, colour = "black") +
  geom_segment(data = Q, aes(x = R, xend = R + 0.22 + 2.2 * lod / maxlod, y = peak, yend = peak,
                             colour = trait), linewidth = 0.18) +
  geom_point(data = Q, aes(x = R + 0.22 + 2.2 * lod / maxlod, y = peak, colour = trait), size = 0.55) +
  geom_segment(data = Q, aes(x = R + 0.22 + 2.2 * lod / maxlod, xend = R + 0.30 + 2.2 * lod / maxlod,
                             y = peak, yend = lab_y), colour = "grey55", linewidth = 0.1) +
  geom_text(data = Q, aes(x = R + 0.34 + 2.2 * lod / maxlod, y = lab_y,
                          label = paste0("q", tr_short, "-", family, "-", chr)),
            hjust = 0, size = 1.25, colour = "grey20") +
  scale_fill_manual(values = TRCOL, name = NULL) +
  scale_colour_manual(values = TRCOL, guide = "none") +
  facet_wrap(~ chr, ncol = 7, scales = "free_y", drop = FALSE) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(breaks = NULL, limits = c(-0.8, 9.0)) +
  theme_void(base_size = 6) +
  theme(strip.text = element_text(size = 5.5, face = "bold"),
        panel.spacing = unit(2, "pt"), legend.position = "bottom",
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = 5))
ggsave(file.path(FIG, "Fig2c_metaQTL_map.png"), p, width = 7.4, height = 5.2, dpi = 600)
cat("QTL with CI:", nrow(Q), "\n")