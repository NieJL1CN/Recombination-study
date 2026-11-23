suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
OUT <- file.path(ROOT, "Figures/main")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        strip.text = element_text(size = 7),
        plot.margin = margin(3, 4, 3, 3)))

## ---------- Fig 1 ----------
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))
ld <- ld[group %in% c("carrier","non")][order(bin)]
p1a <- ggplot(ld, aes(bin, r2, colour = group)) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.6) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = c(13, 15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) +
  theme(legend.position = c(0.75, 0.75))
reg <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))
regm <- reg[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1b <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  labs(x = NULL, y = "Crossovers per RIL") +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.2)
chr <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
M <- merge(chr, len, by = "chr"); F <- M[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1c <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 500, y = max(F$co)*0.95, label = "r = 0.76", size = 2.4, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))
p1d <- ggplot(viol[group == "landrace"], aes(pop, co)) +
  geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
fig1 <- (p1a | p1b) / (p1c | p1d) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 165, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 165, units = "mm", dpi = 600)
cat("Fig1 done\n")

## ---------- Fig 2 ----------
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv"))
mb <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                 Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8))
eff <- merge(eff, mb, by = "qtl")
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
eff[, chr := factor(chr, levels = chrlev)]
p2a <- ggplot(eff, aes(Mb, chr, colour = direction)) +
  geom_point(aes(size = r2), alpha = 0.85) +
  scale_size_continuous(range = c(0.8, 3), name = expression(R^2)) +
  scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Physical position (Mb)", y = NULL) +
  theme(legend.position = "right", legend.key.size = unit(3,"mm"))
eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) +
  xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 0.6) + coord_flip() +
  scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL)
fig2 <- (p2a | p2b) / (p2c | plot_spacer()) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 150, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 150, units = "mm", dpi = 600)
cat("Fig2 done\n")
