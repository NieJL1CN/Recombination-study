## 17_ld_curve.R -- LD decay curve with smoothed fit and r2 = 0.8 emphasis
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
S <- fread(file.path(TAB, "Table7_LD_decay.csv"))
S <- S[dist <= 60 & is.finite(r2)]
fit <- predict(loess(r2 ~ dist, data = S, span = 0.45), newdata = data.frame(dist = seq(0.5, 60, 0.5)))
sm <- data.table(dist = seq(0.5, 60, 0.5), r2 = fit)
cross <- function(thr) { i <- which(sm$r2 <= thr)[1]; if (is.na(i)) return(NA_real_); sm$dist[i] }
x80 <- cross(0.8); x60 <- cross(0.6)
p <- ggplot() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey45", linewidth = 0.25) +
  geom_point(data = S, aes(x = dist, y = r2), colour = "#2C6E9B", size = 0.5, alpha = 0.35) +
  geom_line(data = sm, aes(x = dist, y = r2), colour = "#1F3B73", linewidth = 0.7) +
  geom_vline(xintercept = x80, linetype = "dotted", colour = "#C0392B", linewidth = 0.4) +
  geom_point(data = data.table(x = x80, y = 0.8), aes(x = x, y = y), colour = "#C0392B", size = 1.6) +
  annotate("text", x = x80 + 1.2, y = 0.87,
           label = sprintf("r2 = 0.8 at %.1f Mb", x80), hjust = 0, size = 2.2, colour = "#C0392B") +
  annotate("segment", x = x80, xend = x80, y = 0.8, yend = 0.86, colour = "#C0392B", linewidth = 0.25) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "physical distance (Mb)", y = expression(r^2))
ggsave(file.path(FIG, "Fig3h_LD_decay.png"), p, width = 4.2, height = 3.2, dpi = 600)
cat("r2=0.8 at", round(x80, 2), "Mb; r2=0.6 at", round(x60, 2), "Mb\n")
