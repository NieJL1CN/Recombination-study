## 145_sim_0.1cM_foldchange.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261124)
N <- 3000; M <- 12000; NGEN <- 12; NREP <- 10
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, N); h2 <- rep(4L, N); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(N, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - 0.6 * (D1 + D2)
    o <- order(-fit)[seq_len(N)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  tt
}
grid <- data.table(mult = c(1, 1.5, 2, 2.5, 3))
out <- list()
for (i in seq_len(nrow(grid))) {
  r <- 0.001 * grid$mult[i]
  for (rep in 1:NREP) {
    tr <- sim(r)
    out[[length(out)+1]] <- data.table(mult = grid$mult[i], rep = rep, gen = 1:NGEN, TTdd = tr)
  }
  cat("done", grid$mult[i], "\n")
}
R <- rbindlist(out)
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table128_sim_0.1cM_foldchange.csv")
S <- R[, .(TTdd = mean(TTdd)), by = .(mult, gen)]
CR <- R[order(gen), .(cross = { i <- which(TTdd >= 0.5); if (length(i)) gen[i[1]] else NA_integer_ }), by = .(mult, rep)]
print(CR[, .(mean_cross = mean(cross, na.rm = TRUE), sd = sd(cross, na.rm = TRUE)), by = mult])
