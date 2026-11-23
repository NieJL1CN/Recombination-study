## 163_dominance_field_selection.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261206)
M <- 12000; NGEN <- 12; NREP <- 30; KMAX <- 3000
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
    "Trec_Drec" = t == 2 & d <= 1)
}
sim <- function(r, case) {
  h1 <- rep(1L, 50); h2 <- rep(4L, 50); tt <- numeric(NGEN); n <- numeric(NGEN)
  for (g in 1:NGEN) {
    K <- length(h1); idx <- sample.int(K, M, TRUE)
    g1 <- make_gamete(h1[idx], h2[idx], r); g2 <- make_gamete(h1[idx], h2[idx], r)
    T1 <- g1 <= 2; T2 <- g2 <= 2; D1 <- g1 == 1 | g1 == 3; D2 <- g2 == 1 | g2 == 3
    t <- T1 + T2; d <- D1 + D2
    tt[g] <- mean(t == 2 & d == 0)
    keep <- keep_rule(t, d, case)
    h1 <- g1[keep]; h2 <- g2[keep]; n[g] <- length(h1)
    if (length(h1) == 0) { tt[g:NGEN] <- NA; break }
    if (length(h1) > KMAX) { sel <- sample.int(length(h1), KMAX); h1 <- h1[sel]; h2 <- h2[sel] }
  }
  tt
}
out <- list()
for (case in c("Tdom_Ddom","Tdom_Drec","Trec_Ddom","Trec_Drec")) {
  for (trt in c("baseline","QTL")) {
    r <- if (trt == "baseline") 0.02 else 0.02 * 1.173
    for (rep in 1:NREP) {
      tt <- sim(r, case)
      out[[length(out)+1]] <- data.table(case = case, treat = trt, rep = rep,
        cross = if (any(tt >= 0.5, na.rm=TRUE)) min(which(tt >= 0.5)) else NA_integer_,
        first = if (any(tt > 0, na.rm=TRUE)) min(which(tt > 0)) else NA_integer_,
        p5 = tt[5], p8 = tt[8])
    }
  }
}
R <- rbindlist(out)
print(R[, .(cross = round(mean(cross, na.rm=TRUE),2), first = round(mean(first, na.rm=TRUE),2),
            p5 = round(mean(p5, na.rm=TRUE),3), p8 = round(mean(p8, na.rm=TRUE),3)), by = .(case, treat)][order(case, treat)])
fwrite(R, "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy/Tables/Table145_dominance_cases.csv")
