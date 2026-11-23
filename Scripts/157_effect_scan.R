## 157_effect_scan.R -- usage: Rscript 157_effect_scan.R <effect_pct> <frac>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
eff <- as.numeric(args[1]); frac <- as.numeric(args[2])
set.seed(20261130 + as.integer(eff * 10) + as.integer(frac * 100))
M <- 12000; NGEN <- 10; NREP <- 30
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, K, pen = 0.6) {
  h1 <- rep(1L, K); h2 <- rep(4L, K); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - pen * (D1 + D2)
    o <- order(-fit)[seq_len(K)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  as.integer(if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_)
}
K <- as.integer(M * frac)
out <- list()
for (d in c(2, 3)) {
  for (trt in c("baseline", "QTL")) {
    r <- if (trt == "baseline") d/100 else d * (1 + eff/100) / 100
    for (rep in 1:NREP) out[[length(out)+1]] <- data.table(effect = eff, frac = frac, dist = d, treat = trt, rep = rep, gen = sim(r, K))
  }
}
fwrite(rbindlist(out), sprintf("D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table138_eff%03d_frac%03d.csv", as.integer(eff), as.integer(frac*100)))
cat("done", eff, frac, "\n")
