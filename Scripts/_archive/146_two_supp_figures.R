## 146_two_supp_figures.R -- (a) 0.1 cM fold-change; (b) all QTL at 2 cM
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
## (a)
A <- fread(file.path(ROOT, "Tables/Table128_sim_0.1cM_foldchange.csv"))
A[, mult_f := factor(sprintf("x%.1f", mult), levels = sprintf("x%.1f", c(1,1.5,2,2.5,3)))]
pa <- ggplot(A[, .(TTdd = mean(TTdd)), by = .(mult_f, gen)], aes(gen, TTdd, colour = mult_f)) +
  geom_line(linewidth = 0.45) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = c("x1.0" = "grey25", "x1.5" = "#E08A2E", "x2.0" = "#2E8B57",
                                 "x2.5" = "#2C6E9B", "x3.0" = "#6A4C93"),
                      name = "Recombination\nincrease") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(8, "pt"), legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_0.1cM_foldchange.png"), pa, width = 3.8, height = 2.4, dpi = 600)
## (b) all families at 2 cM
B <- fread(file.path(ROOT, "Tables/Table127_sim_tight_distances.csv"))[dist == 2]
pcts <- c(ParW049 = 17.3, ParW579 = 13.9, ParW729 = 12.3, RobW145 = 6.6, ParW007 = 4.5, ParW352 = 3.5, ParW145 = 2.3)
B[, lab := fifelse(family == "baseline", "No QTL", paste0(family, " (+", pcts[family], "%)"))]
lv <- c("No QTL", paste0(c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
                        " (+", c(17.3,13.9,12.3,6.6,4.5,3.5,2.3), "%)"))
B[, lab := factor(lab, levels = lv)]
pb <- ggplot(B[, .(TTdd = mean(TTdd)), by = .(lab, gen)], aes(gen, TTdd, colour = lab)) +
  geom_line(aes(linewidth = lab == "No QTL")) + geom_point(size = 0.4) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_linewidth_manual(values = c("TRUE" = 0.5, "FALSE" = 0.25), guide = "none") +
  scale_colour_manual(values = setNames(c("grey20","#C0392B","#2C6E9B","#2E8B57","#6A4C93","#E08A2E","#8C564B","#7F7F7F"), lv), name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_continuous(breaks = 1:6, limits = c(1, 6)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(legend.position = "bottom", legend.key.width = unit(7, "pt"),
        legend.text = element_text(size = 5))
ggsave(file.path(ROOT, "Figures/FigS_fixation_allQTL_2cM.png"), pb, width = 3.8, height = 2.7, dpi = 600)
cat("done\n")