suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
img_panel <- function(path) { img <- readPNG(path); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1a <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.7,0.8))
p1b <- img_panel(file.path(OUT, "panels", "Fig1b_native.png"))
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
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
fig1 <- p1b / (p1a | p1c) / (p1d | p1e) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.45, 1, 1))
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 215, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 215, units = "mm", dpi = 600)
## Fig3 fix legend
hd <- data.table(qtl = c("5D","6A","7D"), without = c(3.01,3.08,3.45), with = c(2.96,3.11,3.38))
hd[, qtl := factor(qtl, levels = c("5D","6A","7D"))]; hdd <- melt(hd, id.vars = "qtl", variable.name = "model", value.name = "LOD")
p3a <- ggplot(hdd, aes(qtl, LOD, fill = model)) + geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", LOD)), position = position_dodge(0.6), vjust = -0.4, size = 2.1) +
  scale_fill_manual(values = c(without = "#B0B0B0", with = "#4C72B0"), labels = c("no HD","+ HD")) +
  labs(x = "Family", y = "Peak LOD", fill = NULL) + coord_cartesian(ylim = c(2.9, 3.55)) +
  theme(legend.position = c(0.25, 0.85), legend.key.size = unit(2.5,"mm"))
tr <- fread(file.path(TAB, "Table134_all_trait_allele_effects.csv"))
nice <- c(PH_M_cm="Plant height", GY_Calc_tha="Grain yield", GW_M_g1000grn="Thousand-grain weight",
          GN_Calc_Grnm2="Grains per m2", GWid_M_mm="Grain width", GLng_M_mm="Grain length")
tr[, trait2 := nice[trait]]; tr[, label := paste0(qtl, "  ", trait2)]
tr[, neglog := -log10(p_adj)]; tr <- tr[order(neglog)]; tr[, label := factor(label, levels = label)]
p3b <- ggplot(tr, aes(neglog, label)) + geom_col(fill = "#4C72B0", width = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  annotate("text", x = -log10(0.05)+0.03, y = 1.2, label = "P = 0.05", hjust = 0, size = 2.0, angle = 90) +
  labs(x = expression(-log[10]~"adjusted P"), y = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
fig3 <- (p3a / p3b) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 89, height = 120, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 89, height = 120, units = "mm", dpi = 600)
cat("done\n")
