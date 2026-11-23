## 24d_co_methods.R -- direct vs window CO from one shared parse (15 Mb window)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
WIN <- 15e6; STEP <- 7.5e6
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
POS <- unique(M[, .(marker, chr, ppos)])
out <- list()
for (f in COH) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), POS, by = "marker")
  mk <- mk[order(match(chr, CHR), ppos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  n <- ncol(G); dm <- numeric(n); wm <- numeric(n)
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch); if (length(idx) < 5) next
    X <- G[idx, , drop = FALSE]; s <- mk[idx]
    for (j in seq_len(n)) {
      v <- X[, j]; v <- v[!is.na(v)]
      if (length(v) >= 2) dm[j] <- dm[j] + sum(v[-1] != v[-length(v)])
    }
    mx <- max(s$ppos); nw <- ceiling((mx - WIN) / STEP) + 1
    if (nw < 2) next
    st <- matrix(NA_real_, nw, n)
    for (i in seq_len(nw)) {
      a <- which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)
      if (!length(a)) next
      m <- colMeans(X[a, , drop = FALSE], na.rm = TRUE)
      st[i, ] <- ifelse(m <= 0.05, 0, ifelse(m >= 0.95, 1, NA))
    }
    for (j in seq_len(n)) {
      v <- st[, j]; v <- v[!is.na(v)]
      if (length(v) >= 2) wm[j] <- wm[j] + sum(v[-1] != v[-length(v)])
    }
  }
  out[[f]] <- data.table(dataset = f, id = colnames(G), co_marker = dm, co_window = wm)
  cat("done", f, "\n")
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table17_CO_direct_vs_window_15Mb.csv"))
cat("cohort mean: direct", round(mean(R$co_marker), 1), " window", round(mean(R$co_window), 1), "\n")
L <- melt(R[, .(dataset, id, marker = co_marker, window = co_window)],
          id.vars = c("dataset", "id"), variable.name = "method", value.name = "co")
L[, method := factor(method, levels = c("marker", "window"),
                     labels = c("direct marker count", "sliding window (15 Mb)"))]
p <- ggplot(L, aes(x = method, y = co, fill = method)) +
  geom_violin(colour = NA, alpha = 0.5, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.25) +
  geom_jitter(width = 0.05, size = 0.12, alpha = 0.25, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1.1, colour = "#C0392B") +
  scale_fill_manual(values = c("direct marker count" = "#2C6E9B", "sliding window (15 Mb)" = "#2E8B57"), guide = "none") +
  theme_classic(base_size = 7) + theme(axis.text.x = element_text(size = 6)) +
  labs(x = NULL, y = "CO per individual")
ggsave(file.path(FIG, "Fig3a_CO_method_bar.png"), p, width = 3.2, height = 3.2, dpi = 600)
