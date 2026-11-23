## 166_effect_gain_scan.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261209)
M <- 12000; NGEN <- 12; NREP <- 30; K <- 3000
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, K); h2 <- rep(4L, K); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - 0.6 * (D1 + D2)
    o <- order(-fit)[seq_len(K)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  tt
}
out <- list()
for (eff in c(0, 10, 20, 30, 50, 100)) {
  r <- 0.02 * (1 + eff/100)
  for (rep in 1:NREP) {
    tt <- sim(r)
    out[[length(out)+1]] <- data.table(eff = eff, rep = rep, p2 = tt[2], p3 = tt[3], p4 = tt[4],
      g50 = if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_,
      g99 = if (any(tt >= 0.99)) min(which(tt >= 0.99)) else NA_integer_)
  }
}
R <- rbindlist(out)
S <- R[, .(p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), g50 = mean(g50), g99 = mean(g99)), by = eff][order(eff)]
S[, fold_p2 := p2 / p2[eff == 0]]
print(S)
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table148_effect_gain.csv")
