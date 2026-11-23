## 162_culling_population_scan.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261205)
NGEN <- 10; NREP <- 40
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, M, KMAX) {
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
for (M in c(200, 500, 2000)) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
    for (rep in 1:NREP) {
      tt <- sim(r, M, max(50, M))
      out[[length(out)+1]] <- data.table(M = M, treat = trt, rep = rep,
        cross = if (any(tt >= 0.5, na.rm=TRUE)) min(which(tt >= 0.5)) else NA_integer_,
        p3 = tt[3], p4 = tt[4], p5 = tt[5], p6 = tt[6])
    }
  }
}
R <- rbindlist(out)
print(R[, .(cross = round(mean(cross, na.rm=TRUE),2), se = round(sd(cross, na.rm=TRUE)/sqrt(sum(!is.na(cross))),2),
            p3 = round(mean(p3, na.rm=TRUE),3), p4 = round(mean(p4, na.rm=TRUE),3),
            p5 = round(mean(p5, na.rm=TRUE),3), p6 = round(mean(p6, na.rm=TRUE),3)), by = .(M, treat)][order(M, treat)])
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table144_culling_pop.csv")
