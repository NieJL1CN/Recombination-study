## 160_marker_culling_scheme.R -- breeder-style culling: keep T_ and not DD
suppressPackageStartupMessages(library(data.table))
set.seed(20261203)
M <- 12000; NGEN <- 10; NREP <- 30; KMAX <- 3000
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, 50); h2 <- rep(4L, 50)   # start from TD/td heterozygotes
  tt <- numeric(NGEN); nn <- numeric(NGEN)
  for (g in 1:NGEN) {
    K <- length(h1)
    idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    Tdose <- T1 + T2; Ddose <- D1 + D2
    keep <- Tdose >= 1 & Ddose <= 1          # cull tt and DD
    tt[g] <- mean(Tdose == 2 & Ddose == 0)   # TTdd among all offspring
    h1 <- g1[keep]; h2 <- g2[keep]
    nn[g] <- length(h1)
    if (length(h1) > KMAX) { sel <- sample.int(length(h1), KMAX); h1 <- h1[sel]; h2 <- h2[sel] }
  }
  list(tt = tt, n = nn)
}
out <- list()
for (trt in c("baseline","QTL")) {
  r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
  for (rep in 1:NREP) {
    s <- sim(r)
    out[[length(out)+1]] <- data.table(treat = trt, rep = rep, gen = 1:NGEN, p = s$tt, n = s$n)
  }
}
R <- rbindlist(out)
S <- R[, .(p = mean(p), n = mean(n)), by = .(treat, gen)]
print(S[gen <= 8][order(treat, gen)])
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table142_culling_scheme.csv")
