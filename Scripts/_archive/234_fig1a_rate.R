suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6), legend.title = element_text(size = 6),
  plot.tag = element_text(size = 9, face = "bold"), plot.margin = margin(2,3,2,2)))
R <- fread(file.path(TAB, "Table22b_CO_10Mb_rate.csv"))
R[, binMb := bin/1e6]
A <- R[, .(co = mean(co) * 100), by = .(chr, binMb)]
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
A[, chr := factor(chr, levels = chrlev)]
# region colour by position using reference boundaries
B <- data.table(chr = chrlev,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84),
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200),
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373),
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552))
A <- merge(A, B, by = "chr")
A[, region := fifelse(binMb < b1, "R1",
              fifelse(binMb < b2, "R2a",
              fifelse(binMb < b3, "C",
              fifelse(binMb < b4, "R2b", "R3"))))]
A[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
p1a <- ggplot(A, aes(binMb, co, fill = region)) + geom_col(width = 10) +
  facet_grid(chr ~ ., scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c(C = "#B0B0B0", R1 = "#2E6FA3", R2a = "#6A5199", R2b = "#2F8F5B", R3 = "#E08214")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_y_continuous(breaks = c(0, 10), expand = expansion(mult = c(0, 0))) + scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Physical position (Mb)", y = "Crossovers per 100 RILs per 10 Mb") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(size = 5, angle = 0),
        strip.background = element_blank(), axis.text.x = element_text(size = 4.0), axis.text.y = element_text(size = 3.2),
        axis.title = element_text(size = 6.2), panel.spacing = unit(0.2, "mm"),
        legend.position = "bottom", legend.key.size = unit(2,"mm"))
ggsave(file.path(OUT, "panels", "Fig1a_rate.png"), p1a, width = 52, height = 195, units = "mm", dpi = 600)
cat("range", range(A$co), "n", nrow(A), "\n")
print(A[chr == "5B" & binMb >= 230 & binMb <= 260])
