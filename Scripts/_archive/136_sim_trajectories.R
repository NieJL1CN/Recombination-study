## 136_sim_trajectories.R -- TTdd trajectories for baseline vs QTL, by distance
suppressPackageStartupMessages(library(data.table))
set.seed(20261121)
N <- 3000; M <- 12000; NGEN <- 12; NREP <- 10
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, N); h2 <- rep(4L, N); res <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(N, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - 0.6 * (D1 + D2)
    o <- order(-fit)[seq_len(N)]; h1 <- g1[o]; h2 <- g2[o]
    res[g] <- mean((h1 == 2) & (h2 == 2))
  }
  res
}
fam <- data.table(family = c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
                  pct = c(17.3, 13.9, 12.3, 6.6, 4.5, 3.5, 2.3))
dists <- c(1, 2, 5, 10)
out <- list()
for (d in dists) {
  for (rep in 1:NREP) {
    tr <- sim(d / 100)
    out[[length(out)+1]] <- data.table(family = "baseline", pct = 0, dist = d, rep = rep, gen = 1:NGEN, TTdd = tr)
  }
  for (i in seq_len(nrow(fam))) {
    r <- min(0.5, d * (1 + fam$pct[i]/100) / 100)
    for (rep in 1:NREP) {
      tr <- sim(r)
      out[[length(out)+1]] <- data.table(family = fam$family[i], pct = fam$pct[i], dist = d, rep = rep, gen = 1:NGEN, TTdd = tr)
    }
  }
  cat("done", d, "\n")
}
R <- rbindlist(out)
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table123_sim_trajectories.csv")
cat("rows:", nrow(R), "\n")
