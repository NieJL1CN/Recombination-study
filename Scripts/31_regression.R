## 31_regression.R -- CO ~ length and CO ~ PD, whole chromosome and five regions
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
POS <- unique(M[, .(marker, chr, ppos)])
LEN <- merge(POS, B, by = "chr")
LEN[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
      fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
CHL <- LEN[, .(len = (max(ppos) - min(ppos)) / 1e6), by = chr]
RGL <- LEN[, .(len = (max(ppos) - min(ppos)) / 1e6), by = .(chr, region)]
PDCHR <- merge(unique(M[, .(dataset, marker, chr)]), POS, by = c("marker", "chr"))
PDCHR <- merge(PDCHR, B, by = "chr")
PDCHR[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
         fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
PDCHR[, family := sub("_WATDE.*$", "", dataset)]
PDC <- PDCHR[, .(pd = .N), by = .(family, chr)]
PDR <- PDCHR[, .(pd = .N), by = .(family, chr, region)]
EV <- fread(file.path(TAB, "Table18c_CO_events.csv"))
EV[, family := sub("_WATDE.*$", "", dataset)]
NIND <- EV[, .(n_ind = uniqueN(id)), by = dataset]
COC <- EV[, .(k = .N), by = .(dataset, chr)]
COC <- merge(COC, NIND, by = "dataset")
COC[, co := k / n_ind]
COC[, family := sub("_WATDE.*$", "", dataset)]
COC <- COC[, .(co = mean(co)), by = .(family, chr)]
COR <- EV[, .(k = .N), by = .(dataset, chr, region)]
COR <- merge(COR, NIND, by = "dataset")
COR[, co := k / n_ind]
COR[, family := sub("_WATDE.*$", "", dataset)]
COR <- COR[, .(co = mean(co)), by = .(family, chr, region)]
regs <- c("R1","R2a","C","R2b","R3")
fit_lab <- function(d, xv, yv) {
  m <- lm(d[[yv]] ~ d[[xv]]); b <- coef(m)
  p <- summary(m)$coefficients[2, 4]; r2 <- summary(m)$r.squared
  sprintf("y = %.2f + %.4fx\nR2 = %.2f, P = %s", b[1], b[2], r2, format.pval(p, digits = 2))
}
plot_one <- function(d, xv, yv, xlab, ylab) {
  d <- d[is.finite(get(xv)) & is.finite(get(yv))]
  if (nrow(d) < 5) return(NULL)
  mu <- d[, .(x = mean(get(xv)), y = mean(get(yv))), by = chr]
  lab <- fit_lab(d, xv, yv); xr <- range(d[[xv]]); yr <- range(d[[yv]])
  ggplot(d, aes(x = get(xv), y = get(yv))) +
    geom_point(aes(colour = chr), size = 0.55, alpha = 0.45) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey25", linewidth = 0.45, alpha = 0.15) +
    geom_point(data = mu, aes(x = x, y = y), colour = "#C0392B", size = 0.9, inherit.aes = FALSE) +
    geom_text_repel(data = mu, aes(x = x, y = y, label = chr), size = 1.2, colour = "#C0392B",
                    min.segment.length = Inf, seed = 1, inherit.aes = FALSE) +
    annotate("text", x = xr[2], y = yr[2], label = lab, hjust = 1, vjust = 1, size = 1.5, colour = "grey15") +
    scale_colour_manual(values = rep(c("#2C6E9B", "#6A4C93", "#2E8B57"), 7), guide = "none") +
    theme_classic(base_size = 7) + labs(x = xlab, y = ylab)
}
A <- merge(COC, CHL, by = "chr")
p1 <- plot_one(A, "len", "co", "chromosome length (Mb)", "CO per individual")
ggsave(file.path(FIG, "Fig3e_CO_vs_length_chr.png"), p1, width = 3.3, height = 3, dpi = 600)
Ad <- merge(COC, PDC, by = c("family", "chr"))
p2 <- plot_one(Ad, "pd", "co", "parental polymorphism (PD)", "CO per individual")
ggsave(file.path(FIG, "Fig3e_CO_vs_PD_chr.png"), p2, width = 3.3, height = 3, dpi = 600)
for (rg in regs) {
  a <- merge(COR[region == rg], RGL[region == rg], by = c("chr", "region"))
  p <- plot_one(a, "len", "co", paste0(rg, " length (Mb)"), "CO per individual")
  if (!is.null(p)) ggsave(file.path(FIG, paste0("Fig3e_CO_vs_length_", rg, ".png")), p, width = 3.3, height = 3, dpi = 600)
  ad <- merge(a, PDR[region == rg], by = c("family", "chr", "region"))
  p2 <- plot_one(ad, "pd", "co", paste0("parental polymorphism (", rg, ")"), "CO per individual")
  if (!is.null(p2)) ggsave(file.path(FIG, paste0("Fig3e_CO_vs_PD_", rg, ".png")), p2, width = 3.3, height = 3, dpi = 600)
}
cat("regression figures done\n")