suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
## Fig1a: mean crossovers per individual
ev <- fread(file.path(TAB, "Table18c_CO_events.csv"))
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
nlines <- uniqueN(CO$id)
ev[, bin := floor(ppos/5e7)*50]
ev[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
ev[, chr := factor(chr, levels = chrlev)]
cnt <- ev[, .(n = .N), by = .(chr, bin, region)]
cnt[, mean_co := n / nlines]
p1a <- ggplot(cnt, aes(bin, mean_co, fill = region)) + geom_col(width = 49) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x", strip.position = "left") +
  scale_fill_manual(values = c(C = "#B0B0B0", R1 = "#2E6FA3", R2a = "#6A5199", R2b = "#2F8F5B", R3 = "#E08214")) +
  scale_y_continuous(breaks = c(1,2,3), limits = c(0, 3), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Physical position (Mb)", y = "Crossovers per RIL") +
  theme(strip.text.y.left = element_text(size = 6, angle = 0), strip.background = element_blank(),
        axis.text = element_text(size = 4.0), axis.title = element_text(size = 6.5),
        panel.spacing = unit(0.3, "mm"), legend.position = "bottom", legend.key.size = unit(2,"mm"),
        axis.line = element_line(linewidth = 0.3), panel.border = element_blank())
## Fig1 other panels
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1b <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.4) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.75,0.8))
regm <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1c <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.2) + labs(x = NULL, y = "Crossovers per RIL")
chr <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
len <- data.table(chr = chrlev, Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
F <- merge(chr, len, by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 500, y = max(F$co)*0.95, label = "r = 0.76", size = 2.4, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5))
design <- "abb
acd
aee"
fig1 <- p1a + p1b + p1c + p1d + p1e + plot_layout(design = design, widths = c(0.95, 1, 1), heights = c(1,1,1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 195, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 195, units = "mm", dpi = 600)
cat("Fig1 done\n")

## Fig2a: chromosome bodies
q <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                chr = c("3B","1B","4A","7D","4B","7B","6A","5D","6B"),
                Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8),
                trait = c("distal","distal","distal","distal","distal","centromere","whole","whole","whole"))
len[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
q[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
len[, cen := Mb * 0.5]
p2a <- ggplot() +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = 0, ymax = Mb), fill = "grey85", colour = "grey35", linewidth = 0.25) +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = cen - 4, ymax = cen + 4), fill = "white", colour = NA) +
  geom_point(data = q, aes(idx, Mb, colour = trait), size = 1.1) +
  geom_text_repel(data = q, aes(idx, Mb, label = qtl), size = 2.0, direction = "y", nudge_x = 0.45,
                  segment.size = 0.15, min.segment.length = 0) +
  scale_colour_manual(values = c(distal = "#2E6FA3", centromere = "grey30", whole = "#E08214")) +
  scale_x_continuous(breaks = seq_along(chrlev), labels = chrlev, expand = expansion(add = 0.6)) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = NULL, y = "Physical position (Mb)", colour = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5,"mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
ggsave(file.path(OUT, "panels", "Fig2a_native.png"), p2a, width = 150, height = 85, units = "mm", dpi = 600)
## compose Fig2
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv")); eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2a <- { img <- readPNG(file.path(OUT, "panels", "Fig2a_native.png")); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 1.1) + coord_flip() +
  scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL) + theme(legend.position = "none")
fig2 <- p2a / (p2b | p2c) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.2, 1))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 190, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 190, units = "mm", dpi = 600)
cat("Fig2 done\n")
