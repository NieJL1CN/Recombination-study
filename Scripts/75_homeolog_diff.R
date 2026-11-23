## 75_homeolog_diff.R -- amino-acid differences among homeologs
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
ids <- c("TraesPARA_EIv1.0_1336720","TraesPARA_EIv1.0_1262790","TraesPARA_EIv1.0_1437800",
         "TraesPARA_EIv1.0_2580590","TraesPARA_EIv1.0_2301160","TraesPARA_EIv1.0_2423550")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.pep.all.fa.gz"), "r")
res <- list(); cur <- NULL; buf <- character()
repeat {
  ln <- readLines(con, n = 100000, warn = FALSE)
  if (!length(ln)) break
  for (x in ln) {
    if (startsWith(x, ">")) {
      if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
      gid <- sub(".*gene:([^ ]+).*", "\\1", x)
      cur <- if (gid %in% ids) gid else NULL
      buf <- character()
    } else if (!is.null(cur)) buf <- c(buf, trimws(x))
  }
}
if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
close(con)
cat("sequences:", length(res), "\n")
nw <- function(a, b) {
  A <- strsplit(a, "")[[1]]; B <- strsplit(b, "")[[1]]
  n <- length(A); m <- length(B)
  sc <- matrix(0, n + 1, m + 1)
  sc[, 1] <- -1 * (0:n); sc[1, ] <- -1 * (0:m)
  for (i in 1:n) for (j in 1:m) {
    s <- if (A[i] == B[j]) 1 else -1
    sc[i + 1, j + 1] <- max(sc[i, j] + s, sc[i, j + 1] - 1, sc[i + 1, j] - 1)
  }
  i <- n; j <- m; alA <- character(); alB <- character()
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && sc[i + 1, j + 1] == sc[i, j] + if (A[i] == B[j]) 1 else -1) {
      alA <- c(A[i], alA); alB <- c(B[j], alB); i <- i - 1; j <- j - 1
    } else if (i > 0 && sc[i + 1, j + 1] == sc[i, j + 1] - 1) {
      alA <- c(A[i], alA); alB <- c("-", alB); i <- i - 1
    } else {
      alA <- c("-", alA); alB <- c(B[j], alB); j <- j - 1
    }
  }
  list(a = paste(alA, collapse = ""), b = paste(alB, collapse = ""))
}
diffpos <- function(a, b) {
  al <- nw(a, b)
  x <- strsplit(al$a, "")[[1]]; y <- strsplit(al$b, "")[[1]]
  d <- which(x != y)
  data.table(pos = d, ref = x[d], alt = y[d])
}
r4b <- res[["TraesPARA_EIv1.0_1336720"]]; r4a <- res[["TraesPARA_EIv1.0_1262790"]]; r4d <- res[["TraesPARA_EIv1.0_1437800"]]
cat("RECQL1 lengths:", nchar(r4b), nchar(r4a), nchar(r4d), "\n")
d1 <- diffpos(r4a, r4b); d2 <- diffpos(r4d, r4b)
cat("4A vs 4B differences:", nrow(d1), " 4D vs 4B:", nrow(d2), "\n")
setnames(d1, c("pos","aa_4A","aa_4B")); setnames(d2, c("pos","aa_4D","aa_4B"))
fwrite(d1, file.path(TAB, "Table62_RECQL1_4A_vs_4B.csv"))
fwrite(d2, file.path(TAB, "Table62b_RECQL1_4D_vs_4B.csv"))
print(d1); print(d2)
l7d <- res[["TraesPARA_EIv1.0_2580590"]]; l7a <- res[["TraesPARA_EIv1.0_2301160"]]; l7b <- res[["TraesPARA_EIv1.0_2423550"]]
cat("LSD1 lengths:", nchar(l7d), nchar(l7a), nchar(l7b), "\n")
d3 <- diffpos(l7a, l7d); d4 <- diffpos(l7b, l7d)
cat("7A vs 7D:", nrow(d3), " 7B vs 7D:", nrow(d4), "\n")
setnames(d3, c("pos","aa_7A","aa_7D")); setnames(d4, c("pos","aa_7B","aa_7D"))
fwrite(d3, file.path(TAB, "Table63_LSD1_7A_vs_7D.csv"))
fwrite(d4, file.path(TAB, "Table63b_LSD1_7B_vs_7D.csv"))
print(d3); print(d4)
