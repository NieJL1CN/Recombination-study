## 159_trajectory_population.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261202)
M <- 12000; NGEN <- 10; NREP <- 50
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
  tt
}
K <- as.integer(M * 0.25)
out <- list()
for (trt in c("baseline","QTL")) {
  r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
  for (rep in 1:NREP) {
    tt <- sim(r, K)
    out[[length(out)+1]] <- data.table(treat = trt, rep = rep, gen = 1:NGEN, p = tt)
  }
}
R <- rbindlist(out)
S <- R[, .(p = mean(p), se = sd(p)/sqrt(.N)), by = .(treat, gen)]
print(S[gen <= 6][order(treat, gen)])
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table141_trajectory.csv")
