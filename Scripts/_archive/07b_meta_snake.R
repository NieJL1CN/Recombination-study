## 07b_meta_snake.R -- serpentine capsule map (rounded thick segments) with QTL marks
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
TRC <- c("PH_M_cm" = "#2C6E9B", "GY_Calc_tha" = "#2E8B57", "GW_M_g1000grn" = "#E08A2E")
PK <- fread(file.path(TAB, "Table2c_QTL_family_consensus_vs_original.csv"))
PK[, chr := factor(chr, levels = CHR)]
PK[, col := as.integer(chr)]
PK[, row := (col - 1) %/% 7 + 1]
PK[, idx := (col - 1) %% 7 + 1]
PK[, span := 2]
PK[, xpos := idx + (pos / max(pos)) * 1.4 - 0.7, by = chr]
PK[, y := 4.2 - row]

caps <- data.table(chr = CHR, col = seq_along(CHR))
caps[, row := (col - 1) %/% 7 + 1][, idx := (col - 1) %% 7 + 1][, y := 4.2 - row]

p <- ggplot() +
  geom_segment(data = caps, aes(x = idx - 0.40, xend = idx + 0.40, y = y, yend = y),
               linewidth = 7, lineend = "round", colour = "grey88") +
  geom_segment(data = caps, aes(x = idx - 0.40, xend = idx + 0.40, y = y, yend = y),
               linewidth = 9, lineend = "round", colour = "grey40", alpha = 0.0) +
  geom_point(data = PK, aes(x = xpos, y = y, colour = trait, size = lod), alpha = 0.9, stroke = 0) +
  geom_text(data = PK[lod > 6], aes(x = xpos, y = y + 0.34, label = sub("_WATDE.*$", "", dataset)),
            size = 1.3, colour = "grey15") +
  geom_text(data = caps, aes(x = idx, y = y + 0.42, label = chr), size = 2.1, colour = "grey25") +
  scale_colour_manual(values = TRC, name = NULL) +
  scale_size_continuous(range = c(0.5, 3.4), name = "LOD") +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 6) +
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5),
        plot.margin = margin(8, 4, 4, 4))
ggsave(file.path(FIG, "Fig2b_meta_qtl_snake.png"), p, width = 6.6, height = 2.8, dpi = 600)
cat("snake done", nrow(PK), "\n")