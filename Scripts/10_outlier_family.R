## ============================================================================
## 10_outlier_family.R -- why are ParW007 (3D) and ParW049 (7A) aberrant?
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers
M <- M[!is.na(g_rank)]

## 1. per family x chromosome map length (Axiom.1 panel)
L <- M[panel == "AXIOM.1", .(len = max(g_rank) - min(g_rank), n = .N), by = .(dataset, chr)]
L[, med := median(len), by = chr]
L[, ratio := len / med]
fwrite(L, file.path(TAB, "Table4_map_length_by_family.csv"))
cat("families with any chromosome > 3x median:\n")
print(L[ratio > 3][order(-ratio), .(dataset, chr, len = round(len,1), med = round(med,1), ratio = round(ratio,1))][1:10])

## 2. genotype quality of the two outlier families vs the others
fs <- list.files(GEN, pattern = "axiom_rqtl_genotypes\\.csv$", full.names = TRUE)
Q <- rbindlist(lapply(fs, function(f) {
  x <- read.csv(f, check.names = FALSE, colClasses = "character")
  h <- x[1:2, ]; g <- x[-(1:2), ]
  Mx <- as.matrix(g[, -1, drop = FALSE])
  na <- mean(Mx %in% c("H", "NA", "", "-"))
  hh <- mean(Mx == "H"); ab <- mean(Mx %in% c("A", "B"))
  data.table(dataset = sub("_rqtl_genotypes\\.csv$", "", basename(f)),
             n_ind = nrow(Mx), n_marker = ncol(Mx),
             frac_missing = round(na, 4), frac_het = round(hh, 4), frac_called = round(ab, 4))
}))
Q[, family := sub("_.*$", "", dataset)]
fwrite(Q, file.path(TAB, "Table5_genotype_quality.csv"))
cat("\ngenotype quality (outliers vs others):\n")
print(Q[family %in% c("ParW007", "ParW049") | family == "ParW141"][order(dataset)])
print(Q[!family %in% c("ParW007", "ParW049"), .(mean_missing = round(mean(frac_missing),4),
        mean_het = round(mean(frac_het),4))])

## 3. figure: ratio of chromosome length to the family median
p <- ggplot(L, aes(x = chr, y = ratio)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.25) +
  geom_point(aes(colour = dataset %in% c("ParW007_WATDE0002_watseq_axiom",
                                         "ParW049_WATDE0175_watseq_axiom")), size = 1.1, alpha = 0.8) +
  scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "#C0392B"), guide = "none") +
  theme_classic(base_size = 8) + labs(x = NULL, y = "chromosome length / chromosome median")
ggsave(file.path(FIG, "Fig1e_map_length_outliers.png"), p, width = 6, height = 3, dpi = 500)
cat("\nfigure and tables written\n")
