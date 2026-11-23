## 137_plot_fixation_speed.R -- fixation speed curves
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
R <- fread(file.path(ROOT, "Tables/Table123_sim_trajectories.csv"))
S <- R[, .(TTdd = mean(TTdd), se = sd(TTdd)/sqrt(.N)), by = .(family, pct, dist, gen)]
cr <- S[order(gen), .(cross = { i <- which(TTdd >= 0.5); if (length(i)) gen[i[1]] else NA_integer_ }), by = .(family, pct, dist)]
fwrite(cr, file.path(ROOT, "Tables/Table124_sim_crossing50.csv"))
S[, fam_lab := fifelse(family == "baseline", "No QTL",
               paste0(family, " (+", pct, "%)"))]
lv <- c("No QTL", paste0(c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
                        " (+", c(17.3,13.9,12.3,6.6,4.5,3.5,2.3), "%)"))
S[, fam_lab := factor(fam_lab, levels = lv)]
pal <- c("No QTL" = "grey30", "ParW049 (+17.3%)" = "#C0392B", "ParW579 (+13.9%)" = "#2C6E9B",
         "ParW729 (+12.3%)" = "#2E8B57", "RobW145 (+6.6%)" = "#6A4C93",
         "ParW007 (+4.5%)" = "#E08A2E", "ParW352 (+3.5%)" = "#8C564B",
         "ParW145 (+2.3%)" = "#7F7F7F")
p <- ggplot(S, aes(gen, TTdd, colour = fam_lab)) +
  geom_line(aes(linewidth = fam_lab == "No QTL", alpha = fam_lab == "No QTL")) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.25, colour = "grey40") +
  facet_wrap(~ paste0(dist, " cM"), nrow = 1) +
  scale_linewidth_manual(values = c("TRUE" = 0.8, "FALSE" = 0.35), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.75), guide = "none") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), limits = c(0, 1.02)) +
  labs(x = "Generation of selfing with marker-assisted selection",
       y = "Proportion TTdd (target homozygous, drag removed)") +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 5.2), legend.box = "horizontal")
ggsave(file.path(ROOT, "Figures/FigS_fixation_speed.png"), p, width = 7.6, height = 2.8, dpi = 600)
cat("figure done\n")
print(dcast(cr, family ~ dist, value.var = "cross"))
