## 134_test_additivity.R -- test whether two QTL in a family act additively
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
CO <- fread(file.path(ROOT, "Tables/Table22c_CO_rate_individual.csv"))
CO[, distal := R1 + R3]; CO[, whole := R1 + R2a + C + R2b + R3]
readg <- function(ds) {
  x <- read.csv(file.path(GEN, paste0(ds, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]; rownames(g) <- g[[1]]; g
}
analyse <- function(ds, mk1, mk2, ycol, high1, high2) {
  g <- readg(ds); y <- CO[[ycol]][match(rownames(g), CO$id)]
  a <- g[[mk1]]; b <- g[[mk2]]
  ok <- !is.na(y) & a %in% c("A","B") & b %in% c("A","B")
  hi1 <- ifelse(a == high1, 1L, 0L); hi2 <- ifelse(b == high2, 1L, 0L)
  d <- data.table(y = y[ok], h1 = factor(hi1[ok]), h2 = factor(hi2[ok]))
  m <- lm(y ~ h1 * h2, data = d)
  aov_tab <- anova(m)
  print(d[, .(n = .N, mean = round(mean(y), 2), sd = round(sd(y), 2)), by = .(h1, h2)])
  cat("ANOVA p (h1):", signif(aov_tab["h1","Pr(>F)"], 3),
      " p (h2):", signif(aov_tab["h2","Pr(>F)"], 3),
      " p (interaction):", signif(aov_tab["h1:h2","Pr(>F)"], 3), "\n")
  mm <- d[, .(mean = mean(y)), by = .(h1, h2)]
  base <- mm[h1 == "0" & h2 == "0"]$mean
  e1 <- mm[h1 == "1" & h2 == "0"]$mean - base
  e2 <- mm[h1 == "0" & h2 == "1"]$mean - base
  obs <- mm[h1 == "1" & h2 == "1"]$mean - base
  cat(sprintf("base %.2f | effect1 %+.2f | effect2 %+.2f | expected additive %+.2f | observed both %+.2f\n",
              base, e1, e2, e1 + e2, obs))
  invisible(m)
}
cat("== ParW579 (3B + 4B, both Watkins-high)\n")
analyse("ParW579_WATDE0075_dfw_axiom", "AX-95020695.3B", "AX-94987158.4B", "distal", "B", "B")
cat("\n== ParW729 (1B Watkins-high; 4A Paragon-high)\n")
analyse("ParW729_WATDE0100_watseq_KASP", "BS00012029", "BS00021736", "distal", "B", "A")
