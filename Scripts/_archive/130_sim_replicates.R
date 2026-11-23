## 130_sim_replicates.R -- 10 replicates per parameter combination
suppressPackageStartupMessages(library(data.table))
set.seed(20261117)
N <- 3000; M <- 12000; NGEN <- 15; NREP <- 10
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, N); h2 <- rep(4L, N); res <- list()
  for (g in 1:NGEN) {
    idx <- sample.int(N, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - 0.6 * (D1 + D2)
    o <- order(-fit)[seq_len(N)]; h1 <- g1[o]; h2 <- g2[o]
    res[[g]] <- data.table(gen = g, TTdd = mean((h1 == 2) & (h2 == 2)),
                           Dfreq = mean(c(h1 == 1 | h1 == 3, h2 == 1 | h2 == 3)))
  }
  rbindlist(res)
}
grid <- CJ(base_cM = c(1, 2, 5, 10, 20, 50), mult = c(1, 1.02, 1.17, 2, 5, 10))
out <- list()
for (i in seq_len(nrow(grid))) {
  for (rep in 1:NREP) {
    r <- min(0.5, grid$base_cM[i] * grid$mult[i] / 100)
    s <- sim(r)
    out[[length(out)+1]] <- data.table(base_cM = grid$base_cM[i], mult = grid$mult[i], rep = rep,
      gen_TTdd50 = as.integer(if (any(s$TTdd >= 0.5)) min(s$gen[s$TTdd >= 0.5]) else NA_integer_),
      gen_D10 = as.integer(if (any(s$Dfreq <= 0.10)) min(s$gen[s$Dfreq <= 0.10]) else NA_integer_))
  }
  cat("done", grid$base_cM[i], grid$mult[i], "\n")
}
R <- rbindlist(out)
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table117_sim_replicates.csv")
S <- R[, .(gen_TTdd_mean = mean(gen_TTdd50, na.rm = TRUE), gen_TTdd_sd = sd(gen_TTdd50, na.rm = TRUE),
           gen_D_mean = mean(gen_D10, na.rm = TRUE)), by = .(base_cM, mult)]
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table118_sim_replicates_summary.csv")
print(dcast(S, base_cM ~ mult, value.var = "gen_TTdd_mean"))
