## 15b_co_genome_plot.R -- vertical bars per 10 Mb bin, x = physical, y = chromosome
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
REG <- c("R1" = "#2C6E9B", "R2a" = "#6A4C93", "C" = "#BDBDBD", "R2b" = "#2E8B57", "R3" = "#E08A2E")
D <- fread(file.path(TAB, "Table11_CO_10Mb_bins.csv"))
D[, y := as.numeric(factor(chr, levels = rev(CHR)))]
D[, ymin := y - 0.42]
D[, h := co_per_ind / max(co_per_ind) * 0.8]
D[, ymax := ymin + h]
p <- ggplot(D) +
  geom_rect(aes(xmin = bin * 10, xmax = bin * 10 + 10, ymin = ymin, ymax = ymax, fill = region),
            colour = "grey30", linewidth = 0.03) +
  scale_fill_manual(values = REG, name = NULL) +
  scale_y_continuous(breaks = seq_along(CHR), labels = rev(CHR)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic(base_size = 6) +
  theme(axis.text.y = element_text(size = 4.2), axis.ticks = element_line(linewidth = 0.2),
        legend.position = "bottom", legend.key.size = unit(5, "pt"), legend.text = element_text(size = 5)) +
  labs(x = "physical position (Mb)", y = NULL)
ggsave(file.path(FIG, "Fig3d_CO_genome_distribution.png"), p, width = 6.8, height = 4.2, dpi = 600)
cat("done\n")
