## 118_meta_physical_and_figure.R -- physical meta intervals + Fig3i (all QTL, red = >=2 families)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
M <- fread(file.path(TAB, "Table30_CO_metaQTL.csv"))
## physical interval per cluster from family-specific markers
phys <- function(fam, ch, cm) {
  d <- MK[family == fam & chr == ch]
  if (!nrow(d)) return(NA_real_)
  d[which.min(abs(cM - cm))]$ppos
}
M[, fam_first := tstrsplit(families, ";")[[1]]]
M[, Mb_start := mapply(phys, fam_first, chr, start) / 1e6]
M[, Mb_end := mapply(phys, fam_first, chr, end) / 1e6]
M[, Mb_start := pmin(Mb_start, Mb_end)][, Mb_end := pmax(Mb_start, Mb_end)]
fwrite(M[, .(chr, cluster, n_families, n_traits, n_qtl, start, end, Mb_start = round(Mb_start,1),
             Mb_end = round(Mb_end,1), mean_lod, max_lod, multi, families, traits)],
       file.path(TAB, "Table34_CO_metaQTL_physical.csv"))
## figure: all 63 peaks + red marks for multi-family clusters
P <- fread(file.path(TAB, "Table30b_CO_QTL_suggestive_peaks.csv"))
P[, trait := factor(trait, levels = c("distal","middle","centromere","whole"))]
TCOL <- c(distal = "#2C6E9B", middle = "#6A4C93", centromere = "#7F7F7F", whole = "#E08A2E")
CONS <- L$consensus[, .(marker, chr, cM)]
LEN <- CONS[, .(len = max(cM)), by = chr]
P <- merge(P, LEN, by = "chr")
P[, peak := pmin(pos, len)]
setorder(P, chr, peak)
P[, lab_y := peak][, gap := len * 0.06]
for (ch in unique(P$chr)) {
  idx <- which(P$chr == ch); if (length(idx) < 3) next
  g <- P$gap[idx[1]]; y <- P$peak[idx]
  for (k in seq_along(idx)[-1]) if (y[k] - y[k-1] < g) y[k] <- y[k-1] + g
  P$lab_y[idx] <- y
}
B <- data.table(chr = CHR, x = 0)
p <- ggplot() +
  geom_rect(data = merge(LEN, data.table(chr = CHR), by = "chr"),
            aes(xmin = -0.21, xmax = 0.21, ymin = 0, ymax = len), fill = "white", colour = "black", linewidth = 0.25) +
  geom_segment(data = CONS, aes(x = -0.21, xend = 0.21, y = cM, yend = cM), colour = "grey60", linewidth = 0.04, alpha = 0.6) +
  geom_rect(data = P, aes(xmin = -0.21, xmax = 0.21, ymin = peak - 0.5, ymax = peak + 0.5, fill = trait), alpha = 0.9, colour = NA) +
  geom_segment(data = P, aes(x = 0.21, xend = 0.21 + 1.8 * lod / max(lod), y = peak, yend = peak, colour = trait), linewidth = 0.2) +
  geom_text(data = P, aes(x = 0.24 + 1.8 * lod / max(lod), y = lab_y,
                          label = paste0(trait, "-", family, "-", chr)), hjust = 0, size = 1.1, colour = "grey20") +
  geom_point(data = M[multi == TRUE], aes(x = -0.42, y = (start+end)/2), colour = "#C0392B", size = 0.9) +
  geom_text(data = M[multi == TRUE], aes(x = -0.5, y = (start+end)/2, label = paste0("n=", n_families)),
            hjust = 1, size = 1.0, colour = "#C0392B") +
  scale_fill_manual(values = TCOL, name = NULL) +
  scale_colour_manual(values = TCOL, guide = "none") +
  facet_wrap(~ chr, ncol = 7, scales = "free_y", drop = FALSE) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(breaks = NULL, limits = c(-1.1, 8.6)) +
  theme_void(base_size = 6) +
  theme(strip.text = element_text(size = 5.5, face = "bold"), panel.spacing = unit(2, "pt"),
        legend.position = "bottom", legend.key.size = unit(6, "pt"), legend.text = element_text(size = 5))
ggsave(file.path(FIG, "Fig3i_CO_metaQTL.png"), p, width = 7.8, height = 5.6, dpi = 600)
cat("clusters:", nrow(M), " multi:", sum(M$multi), " peaks:", nrow(P), "\n")
