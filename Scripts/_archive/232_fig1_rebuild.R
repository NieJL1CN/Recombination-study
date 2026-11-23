suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6), legend.title = element_text(size = 6),
  plot.tag = element_text(size = 9, face = "bold"), plot.margin = margin(2,3,2,2)))
chrlen <- c("1A"=594,"1B"=689,"1D"=495,"2A"=780,"2B"=801,"2D"=651,"3A"=717,"3B"=830,"3D"=615,
            "4A"=744,"4B"=673,"4D"=509,"5A"=709,"5B"=713,"5D"=566,"6A"=618,"6B"=720,"6D"=473,
            "7A"=736,"7B"=750,"7D"=638)
chrlev <- names(chrlen)
ev <- fread(file.path(TAB, "Table18c_CO_events.csv"))
nl <- uniqueN(ev$id)
ev[, brk := pmin(floor(ppos/1e7)*10, chrlen[chr] - 10)]
ev[, mid := brk + 5]
ev[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
ev[, chr := factor(chr, levels = chrlev)]
cnt <- ev[, .N, by = .(chr, mid, region)]
cnt[, co := 100 * N / nl]
p1a <- ggplot(cnt, aes(mid, co, fill = region)) + geom_col(width = 10) +
  facet_grid(chr ~ ., scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c(C = "#B0B0B0", R1 = "#2E6FA3", R2a = "#6A5199", R2b = "#2F8F5B", R3 = "#E08214")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Physical position (Mb)", y = "Crossovers per 100 RILs per 10 Mb") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(size = 5, angle = 0),
        strip.background = element_blank(), axis.text.x = element_text(size = 4.2), axis.text.y = element_blank(),
        axis.title = element_text(size = 6.2), panel.spacing = unit(0.2, "mm"),
        legend.position = "bottom", legend.key.size = unit(2,"mm"))
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
F <- merge(chr, data.table(chr = chrlev, Mb = as.numeric(chrlen)), by = "chr")[, .(co = mean(co), Mb = mean(Mb)), by = chr]
p1d <- ggplot(F, aes(Mb, co)) + geom_point(size = 1, colour = "#4C72B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "Chromosome length (Mb)", y = "Crossovers per RIL") +
  annotate("text", x = 480, y = max(F$co)*0.98, label = sprintf("r = %.2f", cor(F$Mb, F$co)), size = 2.2, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
fig <- p1a + p1b + p1c + p1d + p1e + plot_layout(design = "abb
acd
aee", widths = c(0.8, 1, 1), heights = c(1, 1, 1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig, width = 183, height = 220, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig, width = 183, height = 220, units = "mm", dpi = 600)
cat("Fig1 rebuilt\n")
