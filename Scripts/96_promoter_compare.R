## 96_promoter_compare.R -- Paragon vs CS RECQL1 3 kb upstream (k-mer comparison)
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
read_fa <- function(f) {
  l <- readLines(f); paste(l[!startsWith(l, ">")], collapse = "")
}
par <- toupper(read_fa(file.path(TAB, "Table91_Paragon_RECQL1_4B_promoter.fa")))
cs  <- toupper(read_fa(file.path(TAB, "Table91b_CS_RECQL1_4B_promoter.fa")))
cat("Paragon upstream length:", nchar(par), " CS sequence length:", nchar(cs), "\n")
## CS sequence includes gene + 3000 upstream; take first 3000 bp as upstream proxy
cs3 <- substr(cs, 1, 3000)
par3 <- substr(par, 1, 3000)
kmers <- function(s, k = 8) {
  n <- nchar(s) - k + 1
  v <- substring(s, 1:n, 1:n + k - 1)
  unique(v[!grepl("[^ACGT]", v)])
}
kp <- kmers(par3); kc <- kmers(cs3)
inter <- length(intersect(kp, kc)); uni <- length(union(kp, kc))
cat("8-mer Jaccard (Paragon vs CS):", round(inter/uni, 3), "\n")
cat("shared 8-mers:", inter, " Paragon:", length(kp), " CS:", length(kc), "\n")
gc <- function(s) (sum(strsplit(s, "")[[1]] %in% c("G","C"))) / nchar(s)
cat("GC Paragon:", round(gc(par3), 3), " CS:", round(gc(cs3), 3), "\n")
## simple shared-prefix check of first 100 bp
cat("\nfirst 100 bp Paragon:", substr(par3, 1, 100), "\n")
cat("first 100 bp CS     :", substr(cs3, 1, 100), "\n")
