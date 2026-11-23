## 93_RECQL1_full_compare.R -- Paragon vs CS RECQL1 for all three homeologs
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); WD <- file.path(ROOT, "Source/WebData")
read_fa <- function(f) {
  l <- readLines(f); ids <- character(); seqs <- list(); cur <- NULL; buf <- character()
  for (x in l) {
    if (startsWith(x, ">")) { if (!is.null(cur)) seqs[[cur]] <- paste(buf, collapse = ""); cur <- sub("^>", "", x); buf <- character() }
    else buf <- c(buf, trimws(x))
  }
  if (!is.null(cur)) seqs[[cur]] <- paste(buf, collapse = "")
  seqs
}
## Paragon proteins direct from the pep file
want <- c(TraesPARA_EIv1.0_1262790 = "PAR_A", TraesPARA_EIv1.0_1336720 = "PAR_B", TraesPARA_EIv1.0_1437800 = "PAR_D")
con <- gzfile(file.path(WD, "Triticum_aestivum_paragon.GCA949126075v1.pep.all.fa.gz"), "r")
res <- list(); cur <- NULL; buf <- character()
repeat {
  ln <- readLines(con, n = 100000, warn = FALSE); if (!length(ln)) break
  for (x in ln) {
    if (startsWith(x, ">")) {
      if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
      gid <- sub(".*gene:([^ ]+).*", "\\1", x); cur <- if (gid %in% names(want)) gid else NULL; buf <- character()
    } else if (!is.null(cur)) buf <- c(buf, trimws(x))
  }
}
if (!is.null(cur) && length(buf)) res[[cur]] <- paste(buf, collapse = "")
close(con)
par <- res
cs <- read_fa(file.path(TAB, "Table85_CS_RECQL1.faa"))
nw <- function(a, b) {
  A <- strsplit(a, "")[[1]]; B <- strsplit(b, "")[[1]]; n <- length(A); m <- length(B)
  sc <- matrix(0, n + 1, m + 1); sc[, 1] <- -1 * (0:n); sc[1, ] <- -1 * (0:m)
  for (i in 1:n) for (j in 1:m) {
    s <- if (A[i] == B[j]) 1 else -1
    sc[i + 1, j + 1] <- max(sc[i, j] + s, sc[i, j + 1] - 1, sc[i + 1, j] - 1)
  }
  i <- n; j <- m; x <- character(); y <- character()
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && sc[i + 1, j + 1] == sc[i, j] + if (A[i] == B[j]) 1 else -1) {
      x <- c(A[i], x); y <- c(B[j], y); i <- i - 1; j <- j - 1
    } else if (i > 0 && sc[i + 1, j + 1] == sc[i, j + 1] - 1) { x <- c(A[i], x); y <- c("-", y); i <- i - 1 }
    else { x <- c("-", x); y <- c(B[j], y); j <- j - 1 }
  }
  list(a = paste(x, collapse = ""), b = paste(y, collapse = ""))
}
out <- list()
for (k in c("A","B","D")) {
  p <- par[[names(want)[match(paste0("PAR_", k), want)]]]
  c0 <- cs[[paste0("CS_", k)]]
  al <- nw(p, c0); x <- strsplit(al$a, "")[[1]]; y <- strsplit(al$b, "")[[1]]
  d <- which(x != y)
  cat("\n== copy", k, "== Paragon", nchar(p), "aa | CS", nchar(c0), "aa | identity", round(mean(x == y) * 100, 1), "% | differences", length(d), "\n")
  if (length(d)) {
    dd <- data.table(copy = k, pos = d, paragon = x[d], CS = y[d])
    print(dd); out[[k]] <- dd
  }
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table86_RECQL1_Paragon_vs_CS.csv"))
cat("\n-- differences inside helicase core (approx 230-600) --\n")
print(R[pos >= 230 & pos <= 600])
