suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        plot.margin = margin(3, 4, 3, 3)))
hd <- data.table(qtl = c("5D","6A","7D"), without = c(3.01,3.08,3.45), with = c(2.96,3.11,3.38))
hd[, qtl := factor(qtl, levels = c("5D","6A","7D"))]
hdd <- melt(hd, id.vars = "qtl", variable.name = "model", value.name = "LOD")
p3a <- ggplot(hdd, aes(qtl, LOD, fill = model)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(without = "#B0B0B0", with = "#4C72B0"), labels = c("no HD","+ HD")) +
  labs(x = "QTL", y = "Peak LOD", fill = NULL) +
  coord_cartesian(ylim = c(0, 3.6)) + theme(legend.position = c(0.22, 0.85))
tr <- fread(file.path(TAB, "Table134_all_trait_allele_effects.csv"))
nice <- c(PH_M_cm="Plant height", GY_Calc_tha="Grain yield", GW_M_g1000grn="Thousand-grain weight",
          GN_Calc_Grnm2="Grains per m2", GWid_M_mm="Grain width", GLng_M_mm="Grain length")
tr[, trait2 := nice[trait]]; tr[, label := paste0(qtl, "  ", trait2)]
tr[, neglog := -log10(p_adj)]; tr <- tr[order(neglog)]; tr[, label := factor(label, levels = label)]
p3b <- ggplot(tr, aes(neglog, label)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  labs(x = expression(-log[10]~"adjusted P"), y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)))
fig3 <- (p3a / p3b) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 89, height = 115, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 89, height = 115, units = "mm", dpi = 600)

## Fig2 layout refinement
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
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.2)
p2c <- ggplot(eff, aes(qtl, ymin = mean_A, ymax = mean_B, colour = direction)) +
  geom_linerange(linewidth = 0.6) + geom_point(aes(y = mean_A), size = 0.4) + geom_point(aes(y = mean_B), size = 0.4) +
  coord_flip() + scale_colour_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = NULL, y = "Crossovers at peak marker", colour = NULL) + theme(legend.position = "none")
fig2 <- (p2a | p2b) / p2c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 150, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 150, units = "mm", dpi = 600)
cat("done\n")
