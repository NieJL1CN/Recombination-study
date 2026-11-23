## 164_full_fixation.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261207)
M <- 12000; NGEN <- 15; NREP <- 30; KMAX <- 3000
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
sim <- function(r, case) {
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
  tt
}
out <- list()
for (case in c("Tdom_Ddom","Tdom_Drec","Trec_Ddom","Trec_Drec","marker_cull")) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
    for (rep in 1:NREP) {
      tt <- sim(r, case)
      g99 <- if (any(tt >= 0.99, na.rm=TRUE)) min(which(tt >= 0.99)) else NA_integer_
      g100 <- if (any(tt >= 1, na.rm=TRUE)) min(which(tt >= 1)) else NA_integer_
      out[[length(out)+1]] <- data.table(case = case, treat = trt, rep = rep, g99 = g99, g100 = g100,
        p10 = tt[10], p15 = tt[15])
    }
  }
}
R <- rbindlist(out)
print(R[, .(g99 = round(mean(g99, na.rm=TRUE),2), g100 = round(mean(g100, na.rm=TRUE),2),
            reached100 = mean(!is.na(g100)), p10 = round(mean(p10, na.rm=TRUE),3)), by = .(case, treat)][order(case, treat)])
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table146_full_fixation.csv")
