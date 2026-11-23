## 158_mapping_function_check.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261201)
M <- 12000; NGEN <- 10; NREP <- 30
make_gamete <- function(h1, h2, r) {
  n <- length(h1); pick <- runif(n) < 0.5
  a <- ifelse(pick, h1, h2); b <- ifelse(pick, h2, h1)
  Ta <- a <= 2; Da <- a == 1 | a == 3; Db <- b == 1 | b == 3
  sw <- runif(n) < r
  1L + (1L - Ta) * 2L + (1L - ifelse(sw, Db, Da))
}
sim <- function(r, K, pen = 0.6) {
  h1 <- rep(1L, K); h2 <- rep(4L, K); tt <- numeric(NGEN)
  for (g in 1:NGEN) {
    idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    fit <- (T1 + T2) - pen * (D1 + D2)
    o <- order(-fit)[seq_len(K)]; h1 <- g1[o]; h2 <- g2[o]
    tt[g] <- mean((h1 == 2L) & (h2 == 2L))
  }
  as.integer(if (any(tt >= 0.5)) min(which(tt >= 0.5)) else NA_integer_)
}
haldane <- function(d_cM) { d <- d_cM/100; (1 - exp(-2*d))/2 }
out <- list()
for (frac in c(0.25, 0.50)) {
  K <- as.integer(M * frac)
  for (map in c("linear", "haldane")) {
    for (eff in c(17.3, 50)) {
      d <- 2
      base <- if (map == "linear") d/100 else haldane(d)
      qtl <- base * (1 + eff/100)
      for (trt in c("baseline","QTL")) {
        r <- if (trt == "baseline") base else qtl
        for (rep in 1:NREP) out[[length(out)+1]] <- data.table(frac=frac, map=map, eff=eff, treat=trt, rep=rep, gen=sim(r, K))
      }
    }
  }
}
R <- rbindlist(out)
S <- R[, .(m = mean(gen, na.rm=TRUE), se = sd(gen, na.rm=TRUE)/sqrt(sum(!is.na(gen)))), by = .(frac, map, eff, treat)]
W <- dcast(S, frac + map + eff ~ treat, value.var = c("m","se"))
W[, saving := m_baseline - m_QTL]
print(W[order(frac, map, eff), .(frac, map, eff, base = round(m_baseline,2), qtl = round(m_QTL,2), saving = round(saving,2))])
fwrite(W, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table140_mapping_function_check.csv")
