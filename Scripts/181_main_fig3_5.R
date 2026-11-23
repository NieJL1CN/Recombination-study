suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        plot.margin = margin(3, 4, 3, 3)))

## ---------- Fig 3 ----------
hd <- data.table(qtl = c("5D","6A","7D"), without = c(3.01,3.08,3.45), with = c(2.96,3.11,3.38))
hd[, qtl := factor(qtl, levels = c("5D","6A","7D"))]
hdd <- melt(hd, id.vars = "qtl", variable.name = "model", value.name = "LOD")
p3a <- ggplot(hdd, aes(qtl, LOD, fill = model)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(without = "#B0B0B0", with = "#4C72B0"), labels = c("without heading date","with heading date")) +
  labs(x = "QTL", y = "Peak LOD", fill = NULL) +
  coord_cartesian(ylim = c(2.8, 3.6)) + theme(legend.position = "top")
tr <- fread(file.path(TAB, "Table134_all_trait_allele_effects.csv"))
tr[, label := paste0(qtl, " ", trait)]
tr[, neglog := -log10(p_adj)]
tr <- tr[order(neglog)]
tr[, label := factor(label, levels = label)]
p3b <- ggplot(tr, aes(neglog, label)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  labs(x = expression(-log[10]~"adjusted P"), y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
fig3 <- (p3a / p3b) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 89, height = 110, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 89, height = 110, units = "mm", dpi = 600)
cat("Fig3 done\n")

## ---------- Fig 4 ----------
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))
p4a <- ggplot(viol, aes(group, co, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.65, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(landrace = "#4C72B0", cultivar = "#DD8452")) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(legend.position = "none")
m4 <- fread(file.path(TAB, "Table157_4A_meta_all.csv"))
m4[, pop := factor(pop, levels = m4[order(effect), pop])]
p4b <- ggplot(m4, aes(effect, pop)) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_point(size = 1) + geom_errorbarh(aes(xmin = effect-1.96*se, xmax = effect+1.96*se), height = 0.2, linewidth = 0.3) +
  labs(x = "Effect on distal crossovers (Paragon - cultivar)", y = NULL)
rep9 <- data.table(fam = c("ParW729","ParW007","ParW352","ParW145","ParW209","ParW141","ParW471","ParW034","ParW468"),
                   peaks = c(9,7,6,5,3,1,1,1,1), match = c(2,0,0,0,0,0,0,0,1))
rep9[, fam := factor(fam, levels = rep9[order(-peaks), fam])]
rep9m <- melt(rep9, id.vars = "fam", variable.name = "type", value.name = "n")
p4c <- ggplot(rep9m, aes(fam, n, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(peaks = "#B0B0B0", match = "#4C72B0"), labels = c("peaks (LOD >= 2)","on a 2017 chromosome")) +
  labs(x = NULL, y = "Number of peaks", fill = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
fig4 <- (p4a | p4b) / p4c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig4.tiff"), fig4, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig4.png"), fig4, width = 183, height = 160, units = "mm", dpi = 600)
cat("Fig4 done\n")

## ---------- Fig 5 ----------
sim <- fread(file.path(TAB, "Table127_sim_tight_distances.csv"))
cross <- sim[, .(gen50 = if (any(TTdd >= 0.5)) min(gen[TTdd >= 0.5]) else NA_integer_), by = .(family, dist, rep)]
cm <- cross[, .(g = mean(gen50, na.rm = TRUE)), by = .(family, dist)]
base <- cm[family == "baseline"]
top <- cm[family %in% c("ParW049","ParW579")]
p5a <- ggplot(cm, aes(dist, g, group = family)) +
  geom_line(data = base, colour = "black", linewidth = 0.6) +
  geom_point(data = base, colour = "black", size = 0.8) +
  geom_line(data = top, aes(colour = family), linewidth = 0.5, linetype = 2) +
  geom_point(data = top, aes(colour = family), size = 0.8) +
  scale_x_continuous(trans = "log10", breaks = c(0.1,0.5,1,2,3)) +
  scale_colour_manual(values = c(ParW049 = "#4C72B0", ParW579 = "#DD8452")) +
  labs(x = "Target-drag distance (cM)", y = "Generations to 50% TTdd", colour = NULL)
sch <- data.table(scheme = c("Index","Marker culling","T dom/D dom","T dom/D rec","T rec/D dom","T rec/D rec"),
                  no = c(3.70,13.68,9.69,13.72,2.00,9.70), qtl = c(3.03,13.42,9.70,13.44,2.00,9.69))
schm <- melt(sch, id.vars = "scheme", variable.name = "type", value.name = "g")
schm[, scheme := factor(scheme, levels = sch$scheme)]
p5b <- ggplot(schm, aes(scheme, g, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(no = "#B0B0B0", qtl = "#4C72B0"), labels = c("no QTL","with QTL")) +
  labs(x = NULL, y = "Generations to 99% TTdd", fill = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
N <- 10^seq(2.7, 4.7, length.out = 60)
pop <- rbindlist(lapply(c(0,17,50,100), function(e) data.table(eff = paste0("+",e,"%"), N = N,
  p = 1 - (1 - (0.02*(1+e/100))^2/4)^N)))
p5c <- ggplot(pop, aes(N, p*100, colour = eff)) + geom_line(linewidth = 0.5) +
  scale_x_log10() + scale_colour_manual(values = c("+0%"="#B0B0B0","+17%"="#4C72B0","+50%"="#DD8452","+100%"="#55A868")) +
  labs(x = "F2 plants", y = "Chance of at least one TTdd (%)", colour = NULL) +
  theme(legend.position = c(0.25, 0.75))
fig5 <- (p5a | p5b) / p5c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig5.tiff"), fig5, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig5.png"), fig5, width = 183, height = 160, units = "mm", dpi = 600)
cat("Fig5 done\n")
