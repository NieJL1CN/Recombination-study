## 40_co_meta.R -- meta-QTL analysis of the four CO phenotypes
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
P <- fread(file.path(TAB, "Table20_CO_QTL_profiles.csv"))
P[, family := sub("_WATDE.*$", "", dataset)]
P[, trait := factor(trait, levels = c("distal", "middle", "centromere", "whole"))]
peak <- P[order(chr, pos), {
  n <- .N
  if (n < 3) .SD[0] else {
    ismax <- c(FALSE, diff(diff(lod) < 0) > 0, FALSE)
    idx <- which(ismax & lod >= 2)
    .SD[idx]
  }
}, by = .(family, trait, chr)]
peak <- peak[, .(family, trait, chr, pos, lod)]
peak <- peak[order(-lod)][, head(.SD, 2), by = .(family, trait, chr)]
setorder(peak, chr, pos)
peak[, cluster := {
  cl <- integer(.N); k <- 1L
  for (i in seq_len(.N)) {
    if (i > 1 && pos[i] - pos[i - 1] > 20) k <- k + 1L
    cl[i] <- k
  }
  cl
}, by = chr]
META <- peak[, .(n_families = uniqueN(family), n_traits = uniqueN(trait), n_qtl = .N,
                 start = round(min(pos), 1), end = round(max(pos), 1),
                 mean_lod = round(mean(lod), 2), max_lod = round(max(lod), 2),
                 traits = paste(sort(unique(as.character(trait))), collapse = ";"),
                 families = paste(sort(unique(family)), collapse = ";")), by = .(chr, cluster)]
META[, span := round(end - start, 1)]
fwrite(META, file.path(TAB, "Table30_CO_metaQTL.csv"))
fwrite(peak, file.path(TAB, "Table30b_CO_QTL_suggestive_peaks.csv"))
cat("suggestive peaks:", nrow(peak), " meta clusters:", nrow(META),
    " with >=2 families:", nrow(META[n_families >= 2]), "\n")
## meta-QTL map: y = chromosome, x = cM
M <- META[n_families >= 2][order(match(chr, CHR))]
M[, chr := factor(chr, levels = rev(CHR))]
PK <- peak[, chr := factor(chr, levels = rev(CHR))]
p <- ggplot() +
  geom_rect(data = M, aes(xmin = start, xmax = end, ymin = as.numeric(chr) - 0.35,
                          ymax = as.numeric(chr) + 0.35),
            fill = "grey85", colour = "grey65", linewidth = 0.2) +
  geom_point(data = PK, aes(x = pos, y = as.numeric(chr), colour = trait,
                            size = lod), alpha = 0.75) +
  geom_text(data = M, aes(x = (start + end) / 2, y = as.numeric(chr) + 0.42,
                          label = paste0("n=", n_families)), size = 1.5, colour = "grey30") +
  scale_colour_manual(values = c(distal = "#2C6E9B", middle = "#6A4C93",
                                 centromere = "#BDBDBD", whole = "#E08A2E"), name = NULL) +
  scale_size_continuous(range = c(0.4, 1.8), name = "LOD") +
  scale_y_continuous(breaks = seq_along(CHR), labels = rev(CHR), expand = expansion(add = 0.6)) +
  theme_classic(base_size = 6) +
  theme(axis.text.y = element_text(size = 4.5), legend.position = "bottom",
        legend.key.size = unit(5, "pt"), legend.text = element_text(size = 5),
        legend.box = "horizontal") +
  labs(x = "consensus position (cM)", y = NULL)
## figure is produced by 41_co_meta_map.R (definitive version); this script only writes tables
cat("figure written\n")