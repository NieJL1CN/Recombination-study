## 132_sim_family_effects.R -- simulation using observed QTL effects per family
suppressPackageStartupMessages(library(data.table))
set.seed(20261118)
N <- 3000; M <- 12000; NGEN <- 15; NREP <- 5
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
fam <- data.table(
  family = c("ParW049","ParW579","ParW729","RobW145","ParW007","ParW352","ParW145"),
  pct = c(17.3, 14.4, 14.0, 6.6, 4.5, 3.5, 2.3))
grid <- CJ(family = fam$family, dist_cM = c(1,2,5,10,20))
grid <- merge(grid, fam, by = "family")
out <- list()
for (i in seq_len(nrow(grid))) {
  r <- min(0.5, grid$dist_cM[i] * (1 + grid$pct[i]/100) / 100)
  for (rep in 1:NREP) {
    s <- sim(r)
    out[[length(out)+1]] <- data.table(family = grid$family[i], pct = grid$pct[i],
      dist_cM = grid$dist_cM[i], rep = rep,
      gen = as.integer(if (any(s$TTdd >= 0.5)) min(s$gen[s$TTdd >= 0.5]) else NA_integer_))
  }
  cat("done", grid$family[i], grid$dist_cM[i], "\n")
}
R <- rbindlist(out)
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table120_sim_family_effects.csv")
S <- R[, .(gen_mean = mean(gen, na.rm = TRUE), gen_sd = sd(gen, na.rm = TRUE)), by = .(family, pct, dist_cM)]
## baseline (no change) per distance
B <- data.table(dist_cM = c(1,2,5,10,20), base = c(4,4,3,3,2))
S <- merge(S, B, by = "dist_cM")
S[, saved := base - gen_mean]
fwrite(S, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table121_sim_family_saved.csv")
print(dcast(S, family ~ dist_cM, value.var = "saved"))
