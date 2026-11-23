## ============================================================================
## 23_co_final.R -- direct (marker switch) vs window CO, aligned correctly
##   window: 20 Mb, step 10 Mb, mean-rounded state, X skipped (reference method)
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
WIN <- as.numeric(Sys.getenv("CO_WIN", "20")) * 1e6
STEP <- as.numeric(Sys.getenv("CO_STEP", "10")) * 1e6
cat("window", WIN/1e6, "Mb step", STEP/1e6, "Mb\n")
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

parse_one <- function(f) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  mk <- mk[order(match(chr, CHR), ppos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]      # align AFTER sorting
  list(G = G, mk = mk)
}

co_one <- function(f) {
  o <- parse_one(f); G <- o$G; mk <- o$mk
  if (!nrow(mk)) return(NULL)
  n <- ncol(G)
  co_mark <- numeric(n); co_win <- numeric(n)
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch); if (length(idx) < 5) next
    X <- G[idx, , drop = FALSE]
    ## direct: genotype switches between consecutive markers
    for (j in seq_len(n)) {
      v <- X[, j]; v <- v[!is.na(v)]
      if (length(v) >= 2) co_mark[j] <- co_mark[j] + sum(v[-1] != v[-length(v)])
    }
    ## window: mean-rounded state per window, X skipped
    s <- mk[idx]; mx <- max(s$ppos); nw <- ceiling((mx - WIN) / STEP) + 1
    if (nw < 2) next
    Mx <- matrix(NA_real_, nw, n)
    for (i in seq_len(nw)) {
      a <- which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)
      if (!length(a)) next
      Mx[i, ] <- colMeans(X[a, , drop = FALSE], na.rm = TRUE)
    }
    ## pure window: mean within 5% of 0 or 1 (reference uses round(mean,1) then
    ## keeps only 0/1 states; implement the same tolerance explicitly)
    st <- ifelse(Mx <= 0.05, 0, ifelse(Mx >= 0.95, 1, NA_real_))
    for (j in seq_len(n)) {
      v <- st[, j]; v <- v[!is.na(v)]
      if (length(v) >= 2) co_win[j] <- co_win[j] + sum(v[-1] != v[-length(v)])
    }
  }
  data.table(dataset = f, id = colnames(G), co_marker = co_mark, co_window = co_win)
}

fs <- unique(M$dataset)
res <- list()
for (f in fs) { r <- try(co_one(f), silent = TRUE)
  if (inherits(r, "try-error") || is.null(r)) { cat("skip", f, "\n"); next }
  res[[f]] <- r; cat("done", f, " markers", nrow(parse_one(f)$mk), " mean win", round(mean(r$co_window),1), "\n") }
R <- rbindlist(res)
fwrite(R, file.path(TAB, paste0("Table17_CO_direct_vs_window_", Sys.getenv("CO_WIN", "20"), "Mb.csv")))
S <- R[, .(n_ind = .N, co_marker = round(mean(co_marker), 1), co_window = round(mean(co_window), 1)), by = dataset]
S[, n_marker := nrow(parse_one(dataset)), by = dataset]
fwrite(S, file.path(TAB, "Table17b_CO_family_summary.csv"))
cat("\noverall mean marker:", round(mean(S$co_marker),1), " window:", round(mean(S$co_window),1), "\n")
print(S[order(-co_window)][1:5])