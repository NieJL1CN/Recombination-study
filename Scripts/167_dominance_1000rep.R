## 167_dominance_1000rep.R -- usage: Rscript 167.R <case> <treat>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
case <- args[1]; trt <- args[2]
set.seed(20261210 + sum(utf8ToInt(case)) + ifelse(trt == "QTL", 7, 0))
M <- 12000; NGEN <- 15; NREP <- 1000; KMAX <- 3000
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
keep_rule <- function(t, d, case) {
  switch(case,
    "Tdom_Ddom" = t >= 1 & d == 0,
    "Tdom_Drec" = t >= 1 & d <= 1,
    "Trec_Ddom" = t == 2 & d == 0,
    "Trec_Drec" = t == 2 & d <= 1,
    "marker_cull" = t >= 1 & d <= 1)
}
sim <- function(r) {
  h1 <- rep(1L, 50); h2 <- rep(4L, 50); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    K <- length(h1); idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    t <- T1 + T2; d <- D1 + D2
    tt[g] <- mean(t == 2 & d == 0)
    keep <- keep_rule(t, d, case)
    h1 <- g1[keep]; h2 <- g2[keep]
    if (length(h1) == 0) { tt[g:NGEN] <- NA; break }
    if (length(h1) > KMAX) { sel <- sample.int(length(h1), KMAX); h1 <- h1[sel]; h2 <- h2[sel] }
  }
  c(g50 = if (any(tt >= 0.5, na.rm=TRUE)) min(which(tt >= 0.5)) else NA_integer_,
    g99 = if (any(tt >= 0.99, na.rm=TRUE)) min(which(tt >= 0.99)) else NA_integer_,
    g100 = if (any(tt >= 1, na.rm=TRUE)) min(which(tt >= 1)) else NA_integer_)
}
r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
R <- t(replicate(NREP, sim(r)))
out <- data.table(case = case, treat = trt,
  g50_mean = mean(R[, "g50"], na.rm = TRUE), g50_se = sd(R[, "g50"], na.rm = TRUE)/sqrt(sum(!is.na(R[, "g50"]))),
  g99_mean = mean(R[, "g99"], na.rm = TRUE), g99_se = sd(R[, "g99"], na.rm = TRUE)/sqrt(sum(!is.na(R[, "g99"]))),
  g100_mean = mean(R[, "g100"], na.rm = TRUE), reached100 = mean(!is.na(R[, "g100"])), n = NREP)
fwrite(out, sprintf("D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table149_%s_%s.csv", case, trt))
cat("done", case, trt, "\n")
