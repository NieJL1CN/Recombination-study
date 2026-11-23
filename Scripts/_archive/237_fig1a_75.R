suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
E <- fread(file.path(TAB, "Table18c_CO_events.csv"))
n <- uniqueN(E$id)
A <- E[, .(co = 100 * .N / n), by = .(chr, bin = ppos/1e6)]
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
B <- data.table(chr = chrlev,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84),
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200),
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373),
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552))
A <- merge(A, B, by = "chr")
A[, region := fifelse(bin < b1, "R1", fifelse(bin < b2, "R2a", fifelse(bin < b3, "C", fifelse(bin < b4, "R2b", "R3"))))]
A[, region := factor(region, levels = c("C","R1","R2a","R2b","R3"))]
A[, chr := factor(chr, levels = rev(chrlev))]
p <- ggplot(A, aes(bin, co, fill = region)) + geom_col(width = 7.5) +
  facet_grid(chr ~ ., scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = c(C = "#BDBDBD", R1 = "#2C6E9B", R2a = "#6A4C93", R2b = "#2E8B57", R3 = "#E08A2E"), name = NULL) +
  scale_y_continuous(breaks = c(0,10), expand = expansion(mult = c(0,0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "physical position (Mb)", y = "Crossovers per 100 RILs per 7.5 Mb") +
  theme_classic(base_size = 6) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(size = 4.5, angle = 0), strip.background = element_blank(),
        axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 3.2),
        axis.title = element_text(size = 5.5), panel.spacing = unit(0.15, "mm"),
        legend.position = "bottom", legend.key.size = unit(5, "pt"), legend.text = element_text(size = 5))
ggsave(file.path(OUT, "panels", "Fig1a_7.5Mb.png"), p, width = 2.0, height = 7.1, dpi = 600)
cat("max", round(max(A$co),1), "\n")
