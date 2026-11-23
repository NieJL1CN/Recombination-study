## ============================================================================
## 09_meta_analysis.R -- meta-QTL: cluster per-family QTL on the consensus map
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
PK <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))
PK <- PK[map == "consensus"]
PK[, family := sub("_WATDE.*$", "", dataset)]
PK[, ci := 2 * sqrt(r2 * (1 - r2) / pmax(lod, 0.1))]   # rough 95% CI width proxy (cM)

meta <- PK[order(match(chr, CHR), pos)]
meta[, grp := floor(pos / 50)]
M <- meta[, .(n_families = uniqueN(family),
              n_qtl = .N,
              families = paste(sort(unique(family)), collapse = ";"),
              traits = paste(sort(unique(trait)), collapse = ";"),
              pos_min = round(min(pos), 1), pos_max = round(max(pos), 1),
              pos_meta = round(weighted.mean(pos, pmax(lod, 0.1)), 1),
              max_lod = round(max(lod), 2),
              mean_r2 = round(mean(r2), 3)), by = .(chr, grp)]
M[, grp := NULL]
setorder(M, -n_families, -max_lod)
fwrite(M, file.path(TAB, "Table2f_metaQTL.csv"))
CH <- PK[, .(n_qtl = .N, n_families = uniqueN(family)), by = chr][order(-n_qtl)]
fwrite(CH, file.path(TAB, "Table2g_QTL_per_chromosome.csv"))
print(CH)
cat("meta-QTL clusters:", nrow(M), "\n")
print(M[1:min(12, .N), .(chr, pos_meta, n_families, max_lod, traits)])

## figure: capsule map with meta-QTL bands
M[, col := match(chr, CHR)]
M[, row := (col - 1) %/% 7 + 1][, idx := (col - 1) %% 7 + 1][, y := 4.2 - row]
M[, x0 := idx - 0.4 + 0.8 * pos_min / 500][, x1 := idx - 0.4 + 0.8 * pos_max / 500]
caps <- data.table(chr = CHR); caps[, col := .I][, row := (col - 1) %/% 7 + 1][, idx := (col - 1) %% 7 + 1][, y := 4.2 - row]
p <- ggplot() +
  geom_segment(data = caps, aes(x = idx - 0.40, xend = idx + 0.40, y = y, yend = y),
               linewidth = 7, lineend = "round", colour = "grey90") +
  geom_segment(data = M, aes(x = x0, xend = x1, y = y, yend = y, colour = n_families),
               linewidth = 2.2, lineend = "round") +
  geom_text(data = caps, aes(x = idx, y = y + 0.42, label = chr), size = 2, colour = "grey25") +
  scale_colour_gradientn(colours = c("#F7E6A3", "#9CCB8B", "#2E8B57", "#1F3B73"),
                         name = "families", breaks = c(1, 2, 3)) +
  coord_cartesian(clip = "off") + theme_void(base_size = 6) +
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),
        plot.margin = margin(8, 4, 4, 4))
ggsave(file.path(FIG, "Fig2c_metaQTL_map.png"), p, width = 6.6, height = 2.6, dpi = 600)
cat("meta figure done\n")