## 127_sim_recombination_gain.R -- recombination and linkage-drag removal under marker-assisted selection
suppressPackageStartupMessages(library(data.table))
set.seed(20261116)
N <- 3000; M <- 12000; NGEN <- 20
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r) {
  h1 <- rep(1L, N); h2 <- rep(4L, N)
  out <- list()
  for (g in 1:NGEN) {
    idx <- sample.int(N, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r)
    g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2
    D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - 0.6 * (D1 + D2)
    o <- order(-fit)[seq_len(N)]
    h1 <- g1[o]; h2 <- g2[o]
    out[[g]] <- data.table(gen = g, TT = mean((h1 <= 2) & (h2 <= 2)),
                           Dfreq = mean(c(h1 == 1 | h1 == 3, h2 == 1 | h2 == 3)),
                           TTdd = mean((h1 == 2) & (h2 == 2)),
                           fit = mean(fit[o]))
  }
  rbindlist(out)[, r := r]
}
out <- rbindlist(lapply(c(0.01, 0.05, 0.10, 0.20, 0.50), sim))
fwrite(out, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table110_sim_recombination_gain.csv")
summ <- out[, .(gen_TTdd50 = as.integer(if (any(TTdd >= 0.5)) min(gen[TTdd >= 0.5]) else NA_integer_),
                gen_Dfreq10 = as.integer(if (any(Dfreq <= 0.10)) min(gen[Dfreq <= 0.10]) else NA_integer_),
                fit_gen5 = round(fit[gen == 5], 3), fit_gen10 = round(fit[gen == 10], 3)), by = r]
print(summ)
