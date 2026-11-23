suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        plot.margin = margin(3, 4, 3, 3)))
sim <- fread(file.path(TAB, "Table127_sim_tight_distances.csv"))
cross <- sim[, .(gen50 = if (any(TTdd >= 0.5)) min(gen[TTdd >= 0.5]) else NA_integer_), by = .(family, dist, rep)]
cm <- cross[, .(g = mean(gen50, na.rm = TRUE)), by = .(family, dist)]
base <- cm[family == "baseline"]
top <- cm[family %in% c("ParW049","ParW579")]
p5a <- ggplot(cm, aes(dist, g)) +
  geom_line(data = base, colour = "black", linewidth = 0.6) +
  geom_point(data = base, colour = "black", size = 0.8) +
  geom_line(data = top, aes(colour = family), linewidth = 0.5, linetype = 2) +
  geom_point(data = top, aes(colour = family), size = 0.8) +
  scale_x_continuous(trans = "log10", breaks = c(0.1,0.5,1,2,3)) +
  scale_colour_manual(values = c(ParW049 = "#4C72B0", ParW579 = "#DD8452")) +
  labs(x = "Target-drag distance (cM)", y = "Generations to 50% TTdd", colour = NULL) +
  theme(legend.position = c(0.8, 0.8))
sch <- data.table(scheme = c("Index","Marker culling","T dom/D dom","T dom/D rec","T rec/D dom","T rec/D rec"),
                  no = c(3.70,13.68,9.69,13.72,2.00,9.70), qtl = c(3.03,13.42,9.70,13.44,2.00,9.69))
schm <- melt(sch, id.vars = "scheme", variable.name = "type", value.name = "g")
schm[, scheme := factor(scheme, levels = sch$scheme)]
p5b <- ggplot(schm, aes(scheme, g, fill = type)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(no = "#B0B0B0", qtl = "#4C72B0"), labels = c("no QTL","with QTL")) +
  labs(x = NULL, y = "Generations to 99% TTdd", fill = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
N <- 10^seq(2.7, 4.7, length.out = 60)
pop <- rbindlist(lapply(c(0,17,50,100), function(e) data.table(eff = factor(paste0("+",e,"%"), levels = c("+0%","+17%","+50%","+100%")), N = N, p = 1-(1-(0.02*(1+e/100))^2/4)^N)))
p5c <- ggplot(pop, aes(N, p*100, colour = eff)) + geom_line(linewidth = 0.5) +
  scale_x_log10() + scale_colour_manual(values = c("+0%"="#B0B0B0","+17%"="#4C72B0","+50%"="#DD8452","+100%"="#55A868")) +
  labs(x = "F2 plants", y = "Chance of at least one TTdd (%)", colour = NULL) +
  theme(legend.position = c(0.25, 0.75))
fig5 <- (p5a | p5b) / p5c + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig5.tiff"), fig5, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig5.png"), fig5, width = 183, height = 160, units = "mm", dpi = 600)
cat("Fig5 done\n")
