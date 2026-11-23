## 170_weight_sensitivity.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261212)
M <- 12000; NGEN <- 10; NREP <- 30; K <- 3000
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, pen) {
  h1 <- rep(1L, K); h2 <- rep(4L, K); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    sc <- (T1 + T2) - pen * (D1 + D2)
    o <- order(-sc)[seq_len(K)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  c(g50 = if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_,
    g99 = if (any(tt >= 0.99)) min(which(tt >= 0.99)) else NA_integer_, p2 = tt[2], p3 = tt[3])
}
out <- list()
for (pen in c(0.55, 0.6, 0.7, 0.8, 0.9)) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
    for (rep in 1:NREP) {
      v <- sim(r, pen)
      out[[length(out)+1]] <- data.table(pen = pen, treat = trt, rep = rep, g50 = v["g50"], g99 = v["g99"], p2 = v["p2"], p3 = v["p3"])
    }
  }
}
R <- rbindlist(out)
S <- R[, .(g50 = mean(g50, na.rm=TRUE), g99 = mean(g99, na.rm=TRUE), p2 = mean(p2), p3 = mean(p3)), by = .(pen, treat)]
print(dcast(S, pen ~ treat, value.var = c("g50","g99","p2","p3")))
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table153_weight_extended.csv")
