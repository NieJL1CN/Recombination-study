## 161_culling_distance_scan.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261204)
M <- 12000; NGEN <- 12; NREP <- 30; KMAX <- 3000
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, 50); h2 <- rep(4L, 50); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    K <- length(h1); idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    Tdose <- T1 + T2; Ddose <- D1 + D2
    keep <- Tdose >= 1 & Ddose <= 1
    tt[g] <- mean(Tdose == 2 & Ddose == 0)
    h1 <- g1[keep]; h2 <- g2[keep]
    if (length(h1) == 0) { tt[g:NGEN] <- NA; break }
    if (length(h1) > KMAX) { sel <- sample.int(length(h1), KMAX); h1 <- h1[sel]; h2 <- h2[sel] }
  }
  tt
}
out <- list()
for (d in c(0.5, 1, 2, 3, 5)) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") d/100 else d * 1.173/100
    for (rep in 1:NREP) {
      tt <- sim(r)
      cross <- if (any(tt >= 0.5, na.rm=TRUE)) min(which(tt >= 0.5)) else NA_integer_
      out[[length(out)+1]] <- data.table(dist = d, treat = trt, rep = rep, cross = cross, p5 = tt[5], p6 = tt[6], p7 = tt[7], p8 = tt[8])
    }
  }
}
R <- rbindlist(out)
print(R[, .(cross = mean(cross, na.rm=TRUE), se = sd(cross, na.rm=TRUE)/sqrt(sum(!is.na(cross))),
            p5 = mean(p5, na.rm=TRUE), p6 = mean(p6, na.rm=TRUE), p7 = mean(p7, na.rm=TRUE), p8 = mean(p8, na.rm=TRUE)), by = .(dist, treat)][order(dist, treat)])
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table143_culling_scan.csv")
