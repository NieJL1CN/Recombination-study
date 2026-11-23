## 92_RECQL1_paragon_vs_CS.R -- compare Paragon and Chinese Spring RECQL1 (A/B/D)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); WD <- file.path(ROOT, "Source/WebData")
read_fa <- function(f) {
  l <- readLines(f)
  ids <- character(); seqs <- character(); cur <- 1L; buf <- character()
  for (x in l) {
    if (startsWith(x, ">")) {
      if (cur <= length(ids)) seqs[cur] <- paste(buf, collapse = "")
      ids <- c(ids, sub("^>", "", x)); buf <- character(); cur <- length(ids) + 1L
    } else buf <- c(buf, trimws(x))
  }
  seqs[length(ids)] <- paste(buf, collapse = "")
  setNames(seqs[seq_along(ids)], ids)
}
par <- read_fa(file.path(TAB, "Table56b_proteins.faa"))
par <- par[c("TraesPARA_EIv1.0_1262790","TraesPARA_EIv1.0_1336720","TraesPARA_EIv1.0_1437800")]
names(par) <- c("PAR_A","PAR_B","PAR_D")
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
res <- list()
for (k in c("A","B","D")) {
  p <- par[[paste0("PAR_", k)]]; c0 <- cs[[paste0("CS_", k)]]
  al <- nw(p, c0)
  x <- strsplit(al$a, "")[[1]]; y <- strsplit(al$b, "")[[1]]
  d <- which(x != y)
  ident <- mean(x == y) * 100
  cat("\n== copy", k, "== Paragon", nchar(p), "aa | CS", nchar(c0), "aa | identity", round(ident, 1), "%\n")
  if (length(d)) {
    dd <- data.table(copy = k, pos = d, paragon = x[d], CS = y[d])
    print(dd)
    res[[k]] <- dd
  } else cat("identical\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table86_RECQL1_Paragon_vs_CS.csv"))
cat("\nhelicase domains (approx, Arabidopsis Q9FT74 215-391 / 413-563 + ~15 aa offset):\n")
if (nrow(R)) print(R[pos >= 230 & pos <= 600])