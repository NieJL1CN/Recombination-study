suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.5, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5),
  plot.margin = margin(2,3,2,2), strip.text = element_text(size = 6.5), strip.background = element_blank()))
## Fig1b native
ev <- fread(file.path(TAB, "Table18c_CO_events.csv"))
ev[, bin := floor(ppos/1e7)*10]
ev[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
cnt <- ev[, .N, by = .(chr, bin, region)]
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
cnt[, chr := factor(chr, levels = chrlev)]
p1b <- ggplot(cnt, aes(bin, N, fill = region)) + geom_col(width = 9) +
  facet_wrap(~ chr, ncol = 7, scales = "free_x") +
  scale_fill_manual(values = c(C = "#B0B0B0", R1 = "#2E6FA3", R2a = "#6A5199", R2b = "#2F8F5B", R3 = "#E08214")) +
  labs(x = "Physical position (Mb)", y = "Crossovers") +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 7), legend.position = "bottom", legend.key.size = unit(2.5,"mm"))
ggsave(file.path(OUT, "panels", "Fig1b_native.png"), p1b, width = 120, height = 95, units = "mm", dpi = 600)
## Fig2a native
len <- data.table(chr = chrlev, Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
q <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                chr = c("3B","1B","4A","7D","4B","7B","6A","5D","6B"),
                Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8),
                trait = c("distal","distal","distal","distal","distal","centromere","whole","whole","whole"))
len[, chr := factor(chr, levels = chrlev)]; q[, chr := factor(chr, levels = chrlev)]
p2a <- ggplot() +
  geom_segment(data = len, aes(x = chr, xend = chr, y = 0, yend = Mb), linewidth = 0.35, colour = "grey45") +
  geom_point(data = q, aes(chr, Mb, colour = trait), size = 1.4) +
  geom_text_repel(data = q, aes(chr, Mb, label = qtl), size = 1.9, direction = "y", nudge_x = 0.35,
                  segment.size = 0.15, min.segment.length = 0) +
  scale_colour_manual(values = c(distal = "#2E6FA3", centromere = "grey40", whole = "#E08214")) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = NULL, y = "Physical position (Mb)", colour = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5,"mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
ggsave(file.path(OUT, "panels", "Fig2a_native.png"), p2a, width = 150, height = 90, units = "mm", dpi = 600)
## Fig3 improved
hd <- data.table(qtl = c("5D","6A","7D"), without = c(3.01,3.08,3.45), with = c(2.96,3.11,3.38))
hd[, qtl := factor(qtl, levels = c("5D","6A","7D"))]; hdd <- melt(hd, id.vars = "qtl", variable.name = "model", value.name = "LOD")
p3a <- ggplot(hdd, aes(qtl, LOD, fill = model)) + geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", LOD)), position = position_dodge(0.65), vjust = -0.35, size = 2.2) +
  scale_fill_manual(values = c(without = "#B0B0B0", with = "#4C72B0"), labels = c("without heading date","with heading date")) +
  labs(x = "Family", y = "Peak LOD", fill = NULL) +
  coord_cartesian(ylim = c(2.9, 3.55)) + theme(legend.position = "top")
tr <- fread(file.path(TAB, "Table134_all_trait_allele_effects.csv"))
nice <- c(PH_M_cm="Plant height", GY_Calc_tha="Grain yield", GW_M_g1000grn="Thousand-grain weight",
          GN_Calc_Grnm2="Grains per m2", GWid_M_mm="Grain width", GLng_M_mm="Grain length")
tr[, trait2 := nice[trait]]; tr[, label := paste0(qtl, "  ", trait2)]
tr[, neglog := -log10(p_adj)]; tr <- tr[order(neglog)]; tr[, label := factor(label, levels = label)]
p3b <- ggplot(tr, aes(neglog, label)) + geom_col(fill = "#4C72B0", width = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  annotate("text", x = -log10(0.05)+0.02, y = 1, label = "P = 0.05", hjust = 0, size = 2.1, angle = 90, vjust = 0) +
  labs(x = expression(-log[10]~"adjusted P"), y = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
fig3 <- (p3a / p3b) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 89, height = 120, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 89, height = 120, units = "mm", dpi = 600)
cat("done\n")
