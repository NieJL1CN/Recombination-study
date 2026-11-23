## 121_meta_sig_physical_fig.R -- Table34 + Fig3i from SIGNIFICANT QTL only
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
S <- fread(file.path(TAB, "Table20b_CO_QTL_peaks.csv"))
S[, family := sub("_WATDE.*$", "", dataset)]
phys <- function(fam, ch, cm) { d <- MK[family == fam & chr == ch]; if (!nrow(d)) return(NA_real_); d[which.min(abs(cM - cm))]$ppos }
M <- fread(file.path(TAB, "Table30_CO_metaQTL.csv"))
M[, fam_first := tstrsplit(families, ";")[[1]]]
M[, Mb := mapply(phys, fam_first, chr, (start+end)/2) / 1e6]
fwrite(M[, .(chr, cluster, n_families, n_traits, n_qtl, start, end, Mb = round(Mb,1),
             mean_lod, max_lod, multi, families, traits)], file.path(TAB, "Table34_CO_metaQTL_physical.csv"))
CONS <- L$consensus[, .(marker, chr, cM)]
LEN <- CONS[, .(len = max(cM)), by = chr]
P <- S[, .(family, trait, chr, pos, lod)]
P <- merge(P, LEN, by = "chr")
P[, peak := pmin(pos, len)]
P[, trait := factor(trait, levels = c("distal","middle","centromere","whole"))]
TCOL <- c(distal = "#2C6E9B", middle = "#6A4C93", centromere = "#7F7F7F", whole = "#E08A2E")
setorder(P, chr, peak)
P[, lab_y := peak][, gap := len * 0.06]
for (ch in unique(P$chr)) {
  idx <- which(P$chr == ch); if (length(idx) < 3) next
  g <- P$gap[idx[1]]; y <- P$peak[idx]
  for (k in seq_along(idx)[-1]) if (y[k] - y[k-1] < g) y[k] <- y[k-1] + g
  P$lab_y[idx] <- y
}
p <- ggplot() +
  geom_rect(data = LEN, aes(xmin = -0.21, xmax = 0.21, ymin = 0, ymax = len), fill = "white", colour = "black", linewidth = 0.25) +
  geom_segment(data = CONS, aes(x = -0.21, xend = 0.21, y = cM, yend = cM), colour = "grey60", linewidth = 0.04, alpha = 0.6) +
  geom_rect(data = P, aes(xmin = -0.21, xmax = 0.21, ymin = peak - 0.8, ymax = peak + 0.8, fill = trait), alpha = 0.95, colour = NA) +
  geom_segment(data = P, aes(x = 0.21, xend = 0.21 + 1.8 * lod / max(lod), y = peak, yend = peak, colour = trait), linewidth = 0.2) +
  geom_text(data = P, aes(x = 0.24 + 1.8 * lod / max(lod), y = lab_y,
                          label = paste0(trait, "-", family, "-", chr)), hjust = 0, size = 1.2, colour = "grey20") +
  scale_fill_manual(values = TCOL, name = NULL) +
  scale_colour_manual(values = TCOL, guide = "none") +
  facet_wrap(~ chr, ncol = 7, scales = "free_y", drop = FALSE) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(breaks = NULL, limits = c(-0.4, 8.6)) +
  theme_void(base_size = 6) +
  theme(strip.text = element_text(size = 5.5, face = "bold"), panel.spacing = unit(2, "pt"),
        legend.position = "bottom", legend.key.size = unit(6, "pt"), legend.text = element_text(size = 5))
ggsave(file.path(FIG, "Fig3i_CO_metaQTL.png"), p, width = 7.6, height = 5.6, dpi = 600)
cat("significant QTL drawn:", nrow(P), " clusters:", nrow(M), "\n")
