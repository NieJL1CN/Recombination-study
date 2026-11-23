## 27_co_heatmap.R -- chromosome x family CO heatmap (complete grid, NA = no markers)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
H <- fread(file.path(TAB, "Table18_CO_by_chr_family.csv"))
PD <- fread(file.path(TAB, "Table19_PD_by_chr.csv"))
G <- CJ(family = unique(H$family), chr = CHR)
G <- merge(G, H, by = c("family", "chr"), all.x = TRUE)
G <- merge(G, PD, by = c("family", "chr"), all.x = TRUE)
G[, chr := factor(chr, levels = rev(CHR))]
ord <- G[, .(m = mean(co, na.rm = TRUE)), by = family][order(-m)]$family
G[, family := factor(family, levels = ord)]
## a chromosome with markers but no state change has CO = 0, not missing
G[is.na(co) & !is.na(pd), co := 0]
G[, status := fifelse(is.na(pd), "no marker", "ok")]
p <- ggplot(G, aes(x = family, y = chr)) +
  geom_tile(data = G[status == "ok"], aes(fill = co), colour = "white", linewidth = 0.3) +
  geom_tile(data = G[status != "ok"], fill = "grey96", colour = "white", linewidth = 0.3) +
  geom_text(data = G[status == "ok"], aes(label = sprintf("%.1f", co)), size = 1.5, colour = "grey20") +
  scale_fill_gradientn(colours = c("#F7F7E8", "#D9EBC0", "#9CCB8B", "#5EA97A", "#2C6E9B"), name = "mean CO") +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 4.5), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
ggsave(file.path(FIG, "Fig3c_CO_heatmap.png"), p, width = 5.4, height = 5.2, dpi = 600)
cat("cells without markers:", G[status == "no marker", .N], " of ", nrow(G), "\n")