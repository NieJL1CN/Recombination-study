## 41_co_meta_map.R -- Fig3i: ONLY significant CO-QTL (LOD > permutation threshold)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
CONS <- L$consensus[, .(marker, chr, cM, panel)]
TCOL <- c(distal = "#2C6E9B", middle = "#6A4C93", centromere = "#7F7F7F", whole = "#E08A2E")
B <- data.table(chr = CHR,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6)
POS <- unique(L$markers[!is.na(g_rank), .(marker, ppos)])
D2 <- merge(CONS, POS, by = "marker")
D2[, midp := (B$b2[match(chr, B$chr)] + B$b3[match(chr, B$chr)]) / 2]
CENT <- D2[, .SD[which.min(abs(ppos - midp[1]))], by = chr][, .(chr, cen = cM)]
LEN <- CONS[, .(len = max(cM)), by = chr]
MKPOS <- CONS[, .(chr, marker, cM)]
## ONLY significant CO-QTL
D46 <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
cat("significant CO-QTL:", nrow(D46), "\n")
Q <- D46[, .(family, trait, chr, pos, lod)]
Q <- merge(Q, LEN, by = "chr")
Q[, trait := factor(trait, levels = c("distal", "middle", "centromere", "whole"))]
Q[, peak := pmin(pos, len)]
## meta clusters among the significant QTL only (>= 2 families within 20 cM)
setorder(Q, chr, trait, peak)
Q[, cluster := {
  cl <- integer(.N); cur <- 1L
  for (i in seq_len(.N)) {
    if (i > 1 && (abs(peak[i] - peak[i - 1]) > 20 || trait[i] != trait[i - 1])) cur <- cur + 1L
    cl[i] <- cur
  }
  cl
}, by = chr]
META <- Q[, .(start = min(peak), end = max(peak), n_families = uniqueN(family),
              trait = trait[1]), by = .(chr, cluster)][n_families >= 2]
cat("meta clusters (>=2 families) among significant QTL:", nrow(META), "\n")
W <- 0.42; R <- W / 2; maxlod <- Q[, max(lod)]
Q[, lab_y := peak][, gap := len * 0.06]
setorder(Q, chr, peak)
for (ch in unique(Q$chr)) {
  idx <- which(Q$chr == ch); if (length(idx) < 2) next
  g <- Q$gap[idx[1]]; y <- Q$peak[idx]
  for (k in seq_along(idx)[-1]) if (y[k] - y[k - 1] < g) y[k] <- y[k - 1] + g
  Q$lab_y[idx] <- y
}
bars <- merge(LEN, CENT, by = "chr")
p <- ggplot() +
  geom_rect(data = bars, aes(xmin = -R, xmax = R, ymin = 0, ymax = len),
            fill = "white", colour = "black", linewidth = 0.25) +
  geom_segment(data = MKPOS, aes(x = -R, xend = R, y = cM, yend = cM),
               colour = "grey55", linewidth = 0.05, alpha = 0.6) +
  geom_point(data = CENT, aes(x = 0, y = cen), size = 0.5, colour = "black") +
  geom_rect(data = Q, aes(xmin = -R, xmax = R, ymin = peak - 0.6, ymax = peak + 0.6, fill = trait),
            alpha = 0.95, colour = NA) +
  geom_segment(data = Q, aes(x = R, xend = R + 0.20 + 2.0 * lod / maxlod, y = peak, yend = peak,
                             colour = trait), linewidth = 0.2) +
  geom_point(data = Q, aes(x = R + 0.20 + 2.0 * lod / maxlod, y = peak, colour = trait), size = 0.65) +
  geom_segment(data = Q, aes(x = R + 0.20 + 2.0 * lod / maxlod, xend = R + 0.28 + 2.0 * lod / maxlod,
                             y = peak, yend = lab_y), colour = "grey55", linewidth = 0.1) +
  geom_text(data = Q, aes(x = R + 0.32 + 2.0 * lod / maxlod, y = lab_y,
                          label = paste0(trait, "-", family, "-", chr)),
            hjust = 0, size = 1.2, colour = "grey20") +
  scale_fill_manual(values = TCOL, name = NULL) +
  scale_colour_manual(values = TCOL, guide = "none") +
  facet_wrap(~ chr, ncol = 7, scales = "free_y", drop = FALSE) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(breaks = NULL, limits = c(-0.9, 9.5)) +
  theme_void(base_size = 6) +
  theme(strip.text = element_text(size = 5.5, face = "bold"),
        panel.spacing = unit(2, "pt"), legend.position = "bottom",
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = 5))
ggsave(file.path(FIG, "Fig3i_CO_metaQTL.png"), p, width = 7.6, height = 5.6, dpi = 600)
cat("Fig3i drawn with", nrow(Q), "significant QTL\n")
