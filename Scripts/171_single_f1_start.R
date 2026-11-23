## 171_single_f1_start.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261213)
M <- 8000; NGEN <- 10; NREP <- 30; FRAC <- 0.25
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, n_f1) {
  h1 <- rep(1L, n_f1); h2 <- rep(4L, n_f1); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    K <- length(h1); idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    sc <- (T1 + T2) - 0.6 * (D1 + D2)
    keep <- as.integer(M * FRAC)
    o <- order(-sc)[seq_len(keep)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  c(g50 = if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_,
    g99 = if (any(tt >= 0.99)) min(which(tt >= 0.99)) else NA_integer_, p2 = tt[2], p3 = tt[3])
}
out <- list()
for (nf1 in c(1, 3000)) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
    for (rep in 1:NREP) {
      v <- sim(r, nf1)
      out[[length(out)+1]] <- data.table(n_f1 = nf1, treat = trt, rep = rep, g50 = v["g50"], g99 = v["g99"], p2 = v["p2"], p3 = v["p3"])
    }
  }
}
R <- rbindlist(out)
print(R[, .(g50 = round(mean(g50, na.rm=TRUE),2), g99 = round(mean(g99, na.rm=TRUE),2),
            p2 = round(mean(p2),4), p3 = round(mean(p3),3)), by = .(n_f1, treat)])
