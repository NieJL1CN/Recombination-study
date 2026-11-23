## ============================================================================
## 06_linkage_heatmap.R -- capsule karyotype with cM/Locus density heat
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
L <- readRDS(file.path(DER, "maps.rds")); CM <- L$consensus
W <- 0.34
cap <- function(cx, cy, w, up = TRUE, n = 60) {
  th <- seq(0, pi, length.out = n)
  data.table(x = cx + w * cos(if (up) rev(th) else th),
             y = cy + (if (up) -1 else 1) * w * sin(if (up) rev(th) else th))
}

for (P in unique(CM$panel)) {
  D <- CM[panel == P, .(chr, marker, cM)]
  if (nrow(D) < 100) next
  D <- D[order(match(chr, CHR), cM)]
  D[, x := match(chr, CHR)]
  D[, gap := shift(cM, type = "lead") - cM, by = chr]
  ends <- D[, .(Lc = max(cM, na.rm = TRUE)), by = chr][, x := match(chr, CHR)]
  G <- D[, {
    Lc <- max(cM, na.rm = TRUE); ys <- seq(0, Lc, by = 1)
    g <- approx(x = cM, y = gap, xout = ys, rule = 2)$y
    g[!is.finite(g)] <- median(gap, na.rm = TRUE)
    g <- frollmean(g, n = 7, align = "center")
    g[!is.finite(g)] <- median(g, na.rm = TRUE)
    .(y = ys, dens = g)
  }, by = chr][, x := match(chr, CHR)]
  capf <- rbind(
    ends[, cbind(cap(x, 0, W, TRUE), dens = G[chr == .BY$chr][1, dens], grp = .BY$chr), by = chr][, chr := NULL],
    ends[, cbind(cap(x, Lc, W, FALSE), dens = G[chr == .BY$chr][.N, dens], grp = .BY$chr), by = chr][, chr := NULL])
  outl <- ends[, {
    top <- cap(x, 0, W, TRUE); bot <- cap(x, Lc, W, FALSE)
    rbind(top, data.table(x = x + W, y = Lc), bot)
  }, by = chr]
  ymax <- max(ends$Lc)
  gm <- median(G$dens, na.rm = TRUE)
  G[!is.finite(dens), dens := gm]
  capf[!is.finite(dens), dens := gm]
  lo <- as.numeric(quantile(G$dens, 0.10, na.rm = TRUE))
  hi <- as.numeric(quantile(G$dens, 0.90, na.rm = TRUE))
  brk <- signif(seq(lo, hi, length.out = 5), 2)
  p <- ggplot() +
    geom_tile(data = G, aes(x = x, y = y, fill = dens), width = 2 * W, height = 1.02) +
    geom_polygon(data = capf, aes(x = x, y = y, fill = dens, group = grp), colour = NA) +
    geom_segment(data = D[!is.na(gap)], aes(x = x - W, xend = x + W, y = cM, yend = cM),
                 colour = "black", linewidth = 0.1) +
    geom_polygon(data = outl, aes(x = x, y = y, group = chr), fill = NA, colour = "black", linewidth = 0.45) +
    scale_fill_gradientn(colours = c("#C0392B", "#E67E22", "#F7E6A3", "#9CCB8B", "#1F3B73"),
                         trans = "log10", limits = c(lo, hi), oob = scales::squish,
                         breaks = brk, labels = brk, name = "Density (cM/Locus)", na.value = "white",
                         guide = guide_colourbar(barwidth = 16, barheight = 0.4, direction = "horizontal",
                                                 title.position = "bottom", ticks.linewidth = 0.25,
                                                 frame.colour = NA)) +
    scale_y_reverse(breaks = seq(0, ymax, by = 10), expand = expansion(mult = c(0.03, 0.03))) +
    scale_x_continuous(position = "top", breaks = seq_along(CHR), labels = CHR,
                       expand = expansion(add = 0.7)) +
    labs(x = NULL, y = "cM") +
    theme_classic(base_size = 7) +
    theme(axis.line.x.top = element_line(linewidth = 0.3),
          axis.line.y.left = element_line(linewidth = 0.3),
          axis.line.x.bottom = element_blank(), axis.line.y.right = element_blank(),
          axis.ticks = element_line(linewidth = 0.25), axis.ticks.length = unit(1.2, "pt"),
          axis.text.x.top = element_text(size = 5.5, angle = 90, vjust = 0, hjust = 0),
          axis.text.y = element_text(size = 5),
          axis.title.y = element_text(size = 6.5),
          panel.grid = element_blank(), panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5),
          legend.margin = margin(1, 1, 1, 1))
  ggsave(file.path(FIG, paste0("Fig1c_linkagemap_", P, ".png")), p, width = 6.6, height = 4.6, dpi = 600)
  cat("done", P, "chr", uniqueN(G$chr), "max cM", ymax, "\n")
}