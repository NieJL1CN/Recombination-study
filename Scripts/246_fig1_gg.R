suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6, colour = "black"),
  axis.title = element_text(size = 7), legend.text = element_text(size = 5.5), legend.title = element_blank(),
  plot.tag = element_text(size = 9, face = "bold"), plot.margin = margin(1,2,1,1)))
CHR <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
A <- fread(file.path(TAB, "Table179_CO_10Mb_marker.csv"))
A[, binMb := bin/1e6]
Bd <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84),
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200),
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373),
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552))
A <- merge(A, Bd, by = "chr")
A[, region := fifelse(binMb < b1, "R1", fifelse(binMb < b2, "R2a", fifelse(binMb < b3, "C", fifelse(binMb < b4, "R2b", "R3"))))]
A[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
A[, chr := factor(chr, levels = CHR)]
p1a <- ggplot(A, aes(binMb, co, fill = region)) +
  geom_col(width = 10, colour = "black", linewidth = 0.06) +
  facet_wrap(~ chr, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = c(C = "#BDBDBD", R1 = "#2C6E9B", R2a = "#6A4C93", R2b = "#2E8B57", R3 = "#E08A2E")) +
  scale_x_continuous(limits = c(0, 840), breaks = c(0, 300, 600)) +
  scale_y_continuous(breaks = c(0, 1, 2), expand = expansion(mult = c(0,0))) +
  labs(x = "Physical position (Mb)", y = "Crossovers per 100 RILs per Mb") +
  theme(strip.text.y.left = element_text(size = 4.4, angle = 0), strip.background = element_blank(),
        panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.25),
        axis.text = element_text(size = 3.6), axis.title = element_text(size = 5.6),
        panel.spacing = unit(0.12, "mm"), legend.position = "bottom", legend.key.size = unit(4.5,"pt"))
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1b <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.3) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.72,0.8))
regm <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1c <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.0) + labs(x = NULL, y = "Crossovers per RIL")
chr <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
F <- merge(chr, len, by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 480, y = max(F$co)*0.98, label = sprintf("r = %.2f", cor(F$Mb, F$co)), size = 2.2, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
fig <- p1a + p1b + p1c + p1d + p1e + plot_layout(design = "abb\nacd\naee", widths = c(1.0, 1, 1), heights = c(1,1,1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig, width = 183, height = 225, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig, width = 183, height = 225, units = "mm", dpi = 600)
cat("done\n")
