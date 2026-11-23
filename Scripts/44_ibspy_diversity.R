## 44_ibspy_diversity.R -- haplotype diversity in CO-QTL intervals (IBSpy)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); DER <- file.path(ROOT, "Source/Derived")
TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
fs <- list.files(WD, pattern = "ibs_haplotypes\\.tsv\\.gz$", full.names = TRUE)
fs <- fs[!grepl("\\(1\\)", fs)]
cat("IBSpy files:", length(fs), "\n")
acc <- list()
for (f in fs) {
  ch <- sub("^chr([0-9][ABD]).*$", "\\1", basename(f))
  d <- fread(f, select = c("chr","start","end","query","group"))
  d[, chr := sub("^chr([0-9][ABD]).*$", "\\1", chr)]
  acc[[ch]] <- d[, .(n_acc = uniqueN(query), n_group = uniqueN(group[group > 0]),
                     shannon = { p <- table(group[group > 0]) / sum(group > 0)
                                 -sum(p * log(p)) }, mean_group = mean(group[group > 0])),
                 by = .(chr, start, end)]
  cat("done", ch, " windows:", nrow(acc[[ch]]), "\n")
}
W <- rbindlist(acc)
fwrite(W, file.path(TAB, "Table33_IBSpy_window_diversity.csv"))
cat("windows total:", nrow(W), "\n")
print(summary(W$shannon))