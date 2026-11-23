## 150_foldchange_breakage.R -- breakage version of the 0.1 cM fold-change figure
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
A <- fread(file.path(ROOT, "Tables/Table128_sim_0.1cM_foldchange.csv"))
A[, mult_f := factor(sprintf("x%.1f", mult), levels = sprintf("x%.1f", c(1,1.5,2,2.5,3)))]
S <- A[, .(cleanT = mean(cleanT, na.rm = TRUE)), by = .(mult_f, gen)]
p <- ggplot(S, aes(gen, cleanT, colour = mult_f)) +
  geom_line(linewidth = 0.45) + geom_point(size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.18, colour = "grey55") +
  scale_colour_manual(values = c("x1.0" = "grey25", "x1.5" = "#E08A2E", "x2.0" = "#2E8B57",
                                 "x2.5" = "#2C6E9B", "x3.0" = "#6A4C93"),
                      name = "Recombination\nincrease") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.02)) +
  labs(x = "Generation", y = "T carriers free of the deleterious allele") +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(8, "pt"), legend.text = element_text(size = 5.5))
ggsave(file.path(ROOT, "Figures/FigS_breakage_0.1cM_foldchange.png"), p, width = 3.8, height = 2.4, dpi = 600)
cat("done\n")
