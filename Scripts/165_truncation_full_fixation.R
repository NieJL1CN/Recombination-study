## 165_truncation_full_fixation.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261208)
M <- 12000; NGEN <- 15; NREP <- 30; K <- 3000
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
for (trt in c("baseline","QTL")) {
  r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
  for (rep in 1:NREP) {
    tt <- sim(r)
    out[[length(out)+1]] <- data.table(treat = trt, rep = rep,
      g99 = if (any(tt >= 0.99, na.rm=TRUE)) min(which(tt >= 0.99)) else NA_integer_,
      g100 = if (any(tt >= 1, na.rm=TRUE)) min(which(tt >= 1)) else NA_integer_,
      p10 = tt[10])
  }
}
R <- rbindlist(out)
print(R[, .(g99 = round(mean(g99, na.rm=TRUE),2), g100 = round(mean(g100, na.rm=TRUE),2),
            reached100 = mean(!is.na(g100)), p10 = round(mean(p10, na.rm=TRUE),3)), by = treat])
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table147_truncation_full.csv")
