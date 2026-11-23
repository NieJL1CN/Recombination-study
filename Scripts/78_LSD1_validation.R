## 78_LSD1_validation.R -- subfamily + developmental association
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
## (a) peak marker allele vs heading date in ParW352
mt <- "AX-94861586-7D"
f <- file.path(GEN, "ParW352_WATDE0047_watseq_axiom_rqtl_genotypes.csv")
x <- fread(f, header = FALSE)
hdr <- as.character(unlist(x[1, ]))
dat <- x[-(1:3)]
setnames(dat, hdr); setnames(dat, 1, "line")
if (mt %in% names(dat)) {
  g <- dat[, c("line", mt), with = FALSE]
  setnames(g, c("line", "gt"))
  g <- g[gt %in% c("A","B")]
  H <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
  m <- merge(g, H[, .(line, Hd = Hd_dto_day)], by = "line")
  m <- m[!is.na(Hd)]
  cat("lines with genotype+Hd:", nrow(m), " A:", sum(m$gt == "A"), " B:", sum(m$gt == "B"), "\n")
  print(m[, .(n = .N, Hd_mean = round(mean(Hd), 2), Hd_sd = round(sd(Hd), 2)), by = gt])
  cat("t-test p:", signif(t.test(Hd ~ gt, data = m)$p.value, 3), "\n")
  cat("Wilcoxon p:", signif(wilcox.test(Hd ~ gt, data = m)$p.value, 3), "\n")
} else cat("marker not found in ParW352 genotype file\n")
## (b) published HD QTL nearest to 393.91 Mb on 7D
P <- fread(file.path(TAB, "Table35_published_QTL_cohort.csv"))
P[, chr2 := sub("_.*$", "", chr)]
H7 <- P[trait %in% c("Ant_dto_day","Hd_dto_day") & chr2 == "7D" & !is.na(bp_peak)]
H7[, dist := round(bp_peak/1e6 - 393.91, 1)]
cat("\npublished HD QTL on 7D (n =", nrow(H7), "), nearest 5:\n")
print(H7[order(abs(dist))][1:5, .(population, trait, bp_Mb = round(bp_peak/1e6,1), dist, LOD = round(LOD,2))])
