suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_pub <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        strip.text = element_text(size = 6), strip.background = element_blank(),
        plot.tag = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.25),
        plot.margin = margin(2,3,2,2))
theme_set(theme_pub)
## Fig4
v <- fread(file.path(TAB, "Table167_co_groups_final.csv")); v[, group := factor(group, levels = c("landrace","cultivar"))]
p4a <- ggplot(v, aes(group, co, fill = group)) + geom_violin(trim = FALSE, alpha = 0.65, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(landrace = "#4C72B0", cultivar = "#DD8452")) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(legend.position = "none") +
  annotate("text", x = 1.5, y = 140, label = "P = 0.0088", size = 2.4)
m4 <- fread(file.path(TAB, "Table157_4A_meta_all.csv"))
m4 <- rbind(m4[, .(pop, effect, se)], data.table(pop = "Pooled", effect = 0.06, se = 0.26))
lv <- c(m4[pop != "Pooled"][order(effect), pop], "Pooled"); m4[, pop := factor(pop, levels = lv)]
p4b <- ggplot(m4, aes(effect, pop)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = effect - 1.96*se, xmax = effect + 1.96*se), orientation = "y", width = 0.18, linewidth = 0.35) +
  geom_point(aes(colour = pop == "Pooled", shape = pop == "Pooled"), size = 1.5) +
  scale_colour_manual(values = setNames(c("black", "#C0392B"), c("FALSE","TRUE")), guide = "none") +
  scale_shape_manual(values = setNames(c(16,18), c("FALSE","TRUE")), guide = "none") +
  labs(x = "Paragon allele effect on distal crossovers", y = NULL) +
  annotate("text", x = 1.1, y = nrow(m4), label = "pooled 0.06 (-0.45 to 0.57)", size = 2.0, hjust = 0)
rep9 <- data.table(fam = c("ParW729","ParW007","ParW352","ParW145","ParW209","ParW141","ParW471","ParW034","ParW468"),
                   peaks = c(9,7,6,5,3,1,1,1,1), match = c(2,0,0,0,0,0,0,0,1))
rep9[, fam := factor(fam, levels = rep9[order(-peaks), fam])]
r9 <- melt(rep9, id.vars = "fam", variable.name = "type", value.name = "n")
p4c <- ggplot(r9, aes(fam, n, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(peaks = "#B0B0B0", match = "#4C72B0"), labels = c("peaks (LOD >= 2)","on a 2017 chromosome")) +
  scale_y_continuous(breaks = 0:9) + annotate("text", x = 9, y = 8.5, label = "3 of 34", size = 2.4, hjust = 0) +
  labs(x = NULL, y = "Number of peaks", fill = NULL) + coord_cartesian(ylim = c(0, 9.2)) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
fig4 <- (p4a | p4b) / p4c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig4.tiff"), fig4, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig4.png"), fig4, width = 183, height = 160, units = "mm", dpi = 600)
## Fig5
sim <- fread(file.path(TAB, "Table127_sim_tight_distances.csv"))
cross <- sim[, .(g50 = if (any(TTdd >= 0.5)) min(gen[TTdd >= 0.5]) else NA_integer_), by = .(family, dist, rep)]
cm <- cross[, .(g = mean(g50, na.rm = TRUE)), by = .(family, dist)]
base <- cm[family == "baseline"]; top <- cm[family %in% c("ParW049","ParW579")]
p5a <- ggplot(cm, aes(dist, g)) +
  geom_line(data = base, colour = "black", linewidth = 0.6) + geom_point(data = base, colour = "black", size = 0.8) +
  geom_line(data = top, aes(colour = family), linewidth = 0.5, linetype = 2) + geom_point(data = top, aes(colour = family), size = 0.8) +
  scale_x_continuous(trans = "log10", breaks = c(0.1,0.5,1,2,3)) +
  scale_colour_manual(values = c(ParW049 = "#4C72B0", ParW579 = "#DD8452")) +
  labs(x = "Target-drag distance (cM)", y = "Generations to 50% TTdd", colour = NULL) + theme(legend.position = c(0.8,0.8))
sch <- data.table(scheme = c("Index","Marker culling","T dom/D dom","T dom/D rec","T rec/D dom","T rec/D rec"),
                  no = c(3.70,13.68,9.69,13.72,2.00,9.70), qtl = c(3.03,13.42,9.70,13.44,2.00,9.69))
sm <- melt(sch, id.vars = "scheme", variable.name = "type", value.name = "g")
sm[, scheme := factor(scheme, levels = sch$scheme)]
p5b <- ggplot(sm, aes(scheme, g, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(no = "#B0B0B0", qtl = "#4C72B0"), labels = c("no QTL","with QTL")) +
  labs(x = NULL, y = "Generations to 99% TTdd", fill = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5))
N <- 10^seq(2.7, 4.7, length.out = 60)
pop <- rbindlist(lapply(c(0,17,50,100), function(e) data.table(eff = factor(paste0("+",e,"%"), levels = c("+0%","+17%","+50%","+100%")), N = N, p = 1-(1-(0.02*(1+e/100))^2/4)^N)))
p5c <- ggplot(pop, aes(N, p*100, colour = eff)) + geom_line(linewidth = 0.5) + scale_x_log10() +
  scale_colour_manual(values = c("+0%"="#B0B0B0","+17%"="#4C72B0","+50%"="#DD8452","+100%"="#55A868")) +
  labs(x = "F2 plants", y = "Chance of at least one TTdd (%)", colour = NULL) + theme(legend.position = c(0.25,0.75))
fig5 <- (p5a | p5b) / p5c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig5.tiff"), fig5, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig5.png"), fig5, width = 183, height = 160, units = "mm", dpi = 600)
cat("Fig4-5 done\n")
