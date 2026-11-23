## 29_co_genome_plot.R -- Fig3d: CO along the genome, 10 Mb bins, region colours
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
REG <- c("R1" = "#2C6E9B", "R2a" = "#6A4C93", "C" = "#BDBDBD", "R2b" = "#2E8B57", "R3" = "#E08A2E")
E0 <- fread(file.path(TAB, "Table18c_CO_events.csv"))
## events whose flanking informative windows are far apart cannot be located
E <- E0[gap_win <= 3]
cat("excluded unlocalised events:", nrow(E0) - nrow(E), "of", nrow(E0), "
")
N <- E[, .(n_ind = uniqueN(id)), by = dataset]
E[, bin := floor(ppos / 10e6)]
C <- E[, .(n = .N), by = .(dataset, chr, bin, region)]
C <- merge(C, N, by = "dataset")
C[, per_ind := n / n_ind]
S <- C[, .(co = mean(per_ind)), by = .(chr, bin, region)][order(match(chr, rev(CHR)), bin)]
S[, y := as.numeric(factor(chr, levels = rev(CHR)))]
S[, h := co / max(co) * 0.8]
p <- ggplot(S) +
  geom_rect(aes(xmin = bin * 10, xmax = bin * 10 + 10, ymin = y - 0.4, ymax = y - 0.4 + h, fill = region),
            colour = "grey35", linewidth = 0.03) +
  scale_fill_manual(values = REG, name = NULL) +
  scale_y_continuous(breaks = seq_along(CHR), labels = rev(CHR)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic(base_size = 6) +
  theme(axis.text.y = element_text(size = 4.2), axis.ticks = element_line(linewidth = 0.2),
        legend.position = "bottom", legend.key.size = unit(5, "pt"), legend.text = element_text(size = 5)) +
  labs(x = "physical position (Mb)", y = NULL)
ggsave(file.path(FIG, "Fig3d_CO_genome_distribution.png"), p, width = 2.6, height = 4.2, dpi = 600)
cat("bins", nrow(S), "\n")