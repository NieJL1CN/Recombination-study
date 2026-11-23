## 151_penalty_sensitivity.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261126)
N <- 3000; M <- 12000; NGEN <- 10; NREP <- 5
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, pen) {
  h1 <- rep(1L, N); h2 <- rep(4L, N); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(N, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - pen * (D1 + D2)
    o <- order(-fit)[seq_len(N)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  as.integer(if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_)
}
out <- list()
for (pen in c(0.55, 0.6, 0.7, 0.8)) {
  for (d in c(1, 2, 3)) {
    for (trt in c("baseline", "ParW049")) {
      r <- if (trt == "baseline") d/100 else d * 1.173 / 100
      for (rep in 1:NREP) out[[length(out)+1]] <- data.table(pen = pen, dist = d, treat = trt, rep = rep, gen = sim(r, pen))
    }
  }
  cat("done", pen, "\n")
}
R <- rbindlist(out)
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table131_penalty_sensitivity.csv")
print(R[, .(mean = mean(gen, na.rm = TRUE), sd = sd(gen, na.rm = TRUE)), by = .(pen, dist, treat)][order(pen, dist, treat)])
