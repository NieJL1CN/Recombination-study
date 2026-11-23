suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
MK <- as.data.table(readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))$markers)
MK <- MK[!is.na(ppos), .(chr, ppos)]
MK[, mb := floor(ppos/5e6)*5]
DEN <- MK[, .(n_mk = .N), by = .(chr, mb)]
E <- fread(file.path(TAB, "Table18c_CO_events.csv"))
E[, n := uniqueN(id), by = dataset]
E[, mb := floor(ppos/5e6)*5]
E <- merge(E, DEN, by = c("chr", "mb"), all.x = TRUE)
E <- E[!is.na(n_mk) & n_mk >= 5]
E[, span := gap_win * 7.5 + 15]
E[, s := pmax(ppos - span/2, 0)][, e := ppos + span/2]
E[, b0 := floor(s/10e6)*10e6][, b1 := floor(e/10e6)*10e6]
B <- E[, .(bin = seq(b0[1], b1[1], by = 10e6)), by = .(dataset, id, chr, s, e, span, n)]
B[, w := pmax(0, pmin(bin + 10e6, e) - pmax(bin, s)) / span / n]
P <- B[, .(co = sum(w)), by = .(dataset, id, chr, bin)][, .(co = mean(co)), by = .(dataset, chr, bin)]
A <- P[, .(co = mean(co) * 100), by = .(chr, bin)]
A[, binMb := bin/1e6]
CHR <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
Bd <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84),
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200),
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373),
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552))
A <- merge(A, Bd, by = "chr")
A[, region := fifelse(binMb < b1, "R1", fifelse(binMb < b2, "R2a", fifelse(binMb < b3, "C", fifelse(binMb < b4, "R2b", "R3"))))]
A[, y := length(CHR) - match(chr, CHR) + 1][, h := co/max(co)*0.8]
p <- ggplot(A) + geom_rect(aes(xmin = binMb, xmax = binMb+10, ymin = y-0.4, ymax = y-0.4+h, fill = region), colour = "grey35", linewidth = 0.02) +
  scale_fill_manual(values = c(C = "#BDBDBD", R1 = "#2C6E9B", R2a = "#6A4C93", R2b = "#2E8B57", R3 = "#E08A2E"), name = NULL) +
  scale_y_continuous(breaks = 1:21, labels = rev(CHR), expand = c(0,1)) +
  theme_classic(base_size = 6) +
  theme(axis.text.y = element_text(size = 4.2), legend.position = "bottom", legend.key.size = unit(5,"pt"), legend.text = element_text(size = 5)) +
  labs(x = "physical position (Mb)", y = NULL)
ggsave(file.path(OUT, "panels", "Fig1a_capsule.png"), p, width = 2.0, height = 7.1, dpi = 600)
cat("max", round(max(A$co),1), "\n")
print(A[chr == "5B" & binMb >= 220 & binMb <= 290, .(binMb, co = round(co,1), region)])
