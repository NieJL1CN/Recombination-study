suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(cowplot)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures"); OUT <- file.path(FIG, "main")
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"), axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9, face = "bold"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7), plot.margin = margin(2, 3, 2, 2)))
img_panel <- function(path) {
  img <- readPNG(path); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void()
}
## Fig1
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1a <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.7,0.8))
p1b <- img_panel(file.path(FIG, "supplementary", "Fig3d_CO_genome_distribution.png"))
regm <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1c <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.2) + labs(x = NULL, y = "Crossovers per RIL")
chr <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
F <- merge(chr, len, by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 500, y = max(F$co)*0.95, label = "r = 0.76", size = 2.4, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]
viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
fig1 <- (p1a | p1b) / (p1c | p1d | p1e) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 180, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 180, units = "mm", dpi = 600)
cat("Fig1 done\n")

## Fig2 with ideogram
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
mb <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                 Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8))
eff <- merge(eff, mb, by = "qtl")
p2a <- img_panel(file.path(OUT, "panels", "metaQTL_ideogram.png"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 0.6) + geom_point(aes(y = mean_A), size = 0.4) + geom_point(aes(y = mean_B), size = 0.4) +
  coord_flip() + scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL) + theme(legend.position = "none")
fig2 <- (p2a | p2b) / p2c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 170, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 170, units = "mm", dpi = 600)
cat("Fig2 done\n")

## Fig4a with P value
viol2 <- fread(file.path(TAB, "Table167_co_groups_final.csv"))
viol2[, group := factor(group, levels = c("landrace","cultivar"))]
p4a <- ggplot(viol2, aes(group, co, fill = group)) + geom_violin(trim = FALSE, alpha = 0.65, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(landrace = "#4C72B0", cultivar = "#DD8452")) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(legend.position = "none") +
  annotate("text", x = 1.5, y = max(viol2$co)*0.95, label = "P = 0.0088", size = 2.6)
m4 <- fread(file.path(TAB, "Table157_4A_meta_all.csv")); m4[, pop := factor(pop, levels = m4[order(effect), pop])]
p4b <- ggplot(m4, aes(effect, pop)) + geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_point(size = 1) + geom_errorbar(aes(xmin = effect-1.96*se, xmax = effect+1.96*se), orientation = "y", width = 0.2, linewidth = 0.3) +
  labs(x = "Effect on distal crossovers (Paragon - cultivar)", y = NULL)
rep9 <- data.table(fam = c("ParW729","ParW007","ParW352","ParW145","ParW209","ParW141","ParW471","ParW034","ParW468"),
                   peaks = c(9,7,6,5,3,1,1,1,1), match = c(2,0,0,0,0,0,0,0,1))
rep9[, fam := factor(fam, levels = rep9[order(-peaks), fam])]
rep9m <- melt(rep9, id.vars = "fam", variable.name = "type", value.name = "n")
p4c <- ggplot(rep9m, aes(fam, n, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(peaks = "#B0B0B0", match = "#4C72B0"), labels = c("peaks (LOD >= 2)","on a 2017 chromosome")) +
  scale_y_continuous(breaks = 0:9) + annotate("text", x = 9, y = 8.5, label = "3 of 34", size = 2.4, hjust = 0) +
  labs(x = NULL, y = "Number of peaks", fill = NULL) + coord_cartesian(ylim = c(0, 9.2)) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
fig4 <- (p4a | p4b) / p4c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig4.tiff"), fig4, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig4.png"), fig4, width = 183, height = 160, units = "mm", dpi = 600)
cat("Fig4 done\n")
