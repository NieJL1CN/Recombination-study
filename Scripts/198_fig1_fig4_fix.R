suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 7, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
## Fig1a: 10 Mb bins, counts / 100
ev <- fread(file.path(TAB, "Table18c_CO_events.csv"))
ev[, bin := floor(ppos/1e7)*10]
ev[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
ev[, chr := factor(chr, levels = chrlev)]
cnt <- ev[, .(n = .N), by = .(chr, bin, region)]
cnt[, co100 := n/100]
p1a <- ggplot(cnt, aes(bin, co100, fill = region)) + geom_col(width = 9) +
  facet_grid(chr ~ ., scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c(C = "#B0B0B0", R1 = "#2E6FA3", R2a = "#6A5199", R2b = "#2F8F5B", R3 = "#E08214")) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
  scale_y_continuous(breaks = c(1,2,3), limits = c(0, 4), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Physical position (Mb)", y = "Crossovers per RIL") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(size = 5.5, angle = 0),
        strip.background = element_blank(), axis.text = element_text(size = 4.6),
        axis.title = element_text(size = 6.5), panel.spacing = unit(0.3, "mm"),
        legend.position = "bottom", legend.key.size = unit(2,"mm"))
## Fig1b-e
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
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"), Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
F <- merge(chr, len, by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
fit <- lm(co ~ Mb, data = F); rr <- cor(F$Mb, F$co)
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1.1, colour = "#4C72B0") +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], linewidth = 0.45, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 480, y = max(F$co)*0.98, label = sprintf("r = %.2f", rr), size = 2.4, hjust = 0)
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
cat("Fig1 done; r =", round(rr,3), "\n")

## Fig4b clearer
m4 <- fread(file.path(TAB, "Table157_4A_meta_all.csv"))
m4 <- rbind(m4[, .(pop, effect, se)], data.table(pop = "Pooled", effect = 0.06, se = 0.26))
lv <- c(m4[pop != "Pooled"][order(effect), pop], "Pooled")
m4[, pop := factor(pop, levels = lv)]
m4[, is_pool := pop == "Pooled"]
p4b <- ggplot(m4, aes(effect, pop)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = effect - 1.96*se, xmax = effect + 1.96*se), orientation = "y", width = 0.18, linewidth = 0.35) +
  geom_point(aes(colour = is_pool, shape = is_pool), size = 1.6) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "#C0392B"), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  labs(x = "Paragon allele effect on distal crossovers", y = NULL) +
  annotate("text", x = 1.1, y = nrow(m4), label = "pooled 0.06 (-0.45 to 0.57)", size = 2.0, hjust = 0)
fig4 <- (function(){
  viol2 <- fread(file.path(TAB, "Table167_co_groups_final.csv")); viol2[, group := factor(group, levels = c("landrace","cultivar"))]
  p4a <- ggplot(viol2, aes(group, co, fill = group)) + geom_violin(trim = FALSE, alpha = 0.65, colour = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(landrace = "#4C72B0", cultivar = "#DD8452")) +
    labs(x = NULL, y = "Crossovers per RIL") + theme(legend.position = "none") +
    annotate("text", x = 1.5, y = 140, label = "P = 0.0088", size = 2.6)
  rep9 <- data.table(fam = c("ParW729","ParW007","ParW352","ParW145","ParW209","ParW141","ParW471","ParW034","ParW468"),
                     peaks = c(9,7,6,5,3,1,1,1,1), match = c(2,0,0,0,0,0,0,0,1))
  rep9[, fam := factor(fam, levels = rep9[order(-peaks), fam])]
  rep9m <- melt(rep9, id.vars = "fam", variable.name = "type", value.name = "n")
  p4c <- ggplot(rep9m, aes(fam, n, fill = type)) + geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(peaks = "#B0B0B0", match = "#4C72B0"), labels = c("peaks (LOD >= 2)","on a 2017 chromosome")) +
    scale_y_continuous(breaks = 0:9) + annotate("text", x = 9, y = 8.5, label = "3 of 34", size = 2.4, hjust = 0) +
    labs(x = NULL, y = "Number of peaks", fill = NULL) + coord_cartesian(ylim = c(0, 9.2)) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  (p4a | p4b) / p4c + plot_annotation(tag_levels = "a")
})()
ggsave(file.path(OUT, "Fig4.tiff"), fig4, width = 183, height = 165, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig4.png"), fig4, width = 183, height = 165, units = "mm", dpi = 600)
cat("Fig4 done\n")
