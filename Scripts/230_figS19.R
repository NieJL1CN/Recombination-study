suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); SUP <- file.path(ROOT, "Figures/supplementary")
M <- fread(file.path(TAB, "Table176_trait_metaQTL.csv"))
CM <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))$consensus
LEN <- CM[, .(len = max(cM)), by = chr]
M <- merge(M, LEN, by = "chr")
M[, tr := c(PH_M_cm = "PH", HD = "HD", Hd_dto_day = "HD", GY_Calc_tha = "GY", GW_M_g1000grn = "GW")[trait]]
M <- M[!is.na(tr)]
CHR <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
M[, col := match(chr, CHR)][, row := (col - 1) %/% 7 + 1][, idx := (col - 1) %% 7 + 1][, y := 4.2 - row]
cap <- data.table(chr = CHR); cap[, col := .I][, row := (col - 1) %/% 7 + 1][, idx := (col - 1) %% 7 + 1][, y := 4.2 - row]
M[, x0 := idx - 0.4 + 0.8 * start / len][, x1 := idx - 0.4 + 0.8 * end / len]
M[, xm := idx - 0.4 + 0.8 * pos / len]
M[, lab := paste0(tr, "-", chr, "-", round(pos))]
M[, k := seq_len(.N), by = chr][, ylab := y - 0.26 - 0.14 * ((k - 1) %% 3)]
p <- ggplot() +
  geom_segment(data = cap, aes(x = idx - 0.40, xend = idx + 0.40, y = y, yend = y), linewidth = 6, lineend = "round", colour = "grey92") +
  geom_segment(data = M, aes(x = x0, xend = pmax(x1, x0 + 0.02), y = y, yend = y, colour = tr), linewidth = 3, lineend = "round") +
  geom_text(data = cap, aes(x = idx, y = y + 0.36, label = chr), size = 2, colour = "grey25") +
  geom_text(data = M, aes(x = xm, y = ylab, label = lab, colour = tr), size = 1.1, angle = 0) +
  scale_colour_manual(values = c(PH = "#2C6E9B", HD = "#6A4C93", GY = "#2E8B57", GW = "#E08A2E"), name = NULL) +
  coord_cartesian(clip = "off") + theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.key.size = unit(5, "pt"),
        plot.margin = margin(6, 4, 4, 4))
ggsave(file.path(SUP, "FigS19_trait_metaQTL_map.png"), p, width = 160, height = 95, units = "mm", dpi = 600)
cat("S19 done, clusters:", nrow(M), "\n")
