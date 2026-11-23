## 26_co_by_chr_region.R -- window CO per chromosome and per region (15 Mb / 7.5 Mb)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN <- file.path(ROOT, "Source/Genotype/Watkins_RIL_Genotypes")
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
WIN <- as.numeric(Sys.getenv("CO_WIN", "5")) * 1e6; STEP <- WIN / 2
SFX <- Sys.getenv("SUFFIX", "")
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, chr, ppos)])

out <- list()
for (f in COH) {
  x <- read.csv(file.path(GEN, paste0(f, "_rqtl_genotypes.csv")), check.names = FALSE, colClasses = "character")
  g <- x[-(1:2), ]
  G <- t(apply(as.matrix(g[, -1, drop = FALSE]), 2, function(z) {
    z <- trimws(z); o <- rep(NA_real_, length(z))
    o[z %in% c("A","a","1","AA")] <- 0; o[z %in% c("B","b","2","BB")] <- 1; o }))
  rownames(G) <- colnames(x)[-1]; colnames(G) <- g[[1]]
  mk <- merge(data.table(marker = rownames(G)), pos, by = "marker")
  mk <- mk[order(match(chr, CHR), ppos)]
  G <- G[match(mk$marker, rownames(G)), , drop = FALSE]
  mk <- merge(mk, B, by = "chr")
  rec <- list()
  for (ch in unique(mk$chr)) {
    idx <- which(mk$chr == ch)
    s <- mk[idx]; if (nrow(s) < 5) next
    mx <- max(s$ppos); nw <- ceiling((mx - WIN) / STEP) + 1
    if (nw < 2) next
    Mx <- matrix(NA_real_, nw, ncol(G)); mid <- numeric(nw)
    for (i in seq_len(nw)) {
      a <- which(s$ppos > (i - 1) * STEP & s$ppos <= (i - 1) * STEP + WIN)
      mid[i] <- (i - 1) * STEP + WIN / 2
      if (!length(a)) next
      Mx[i, ] <- colMeans(G[idx[a], , drop = FALSE], na.rm = TRUE)
    }
    st <- ifelse(Mx <= 0.05, 0, ifelse(Mx >= 0.95, 1, NA_real_))
    reg <- fifelse(mid < s$b1[1], "R1", fifelse(mid < s$b2[1], "R2a",
            fifelse(mid < s$b3[1], "C", fifelse(mid < s$b4[1], "R2b", "R3"))))
    for (j in seq_len(ncol(G))) {
      v <- st[, j]; id2 <- which(!is.na(v))
      if (length(id2) < 2) next
      w <- v[id2]
      k <- which(w[-1] != w[-length(w)])
      if (!length(k)) next
      ## crossover lies between the two informative windows; use the midpoint
      ## so that long X gaps (centromeres) do not push events to one edge
      pos_ev <- (mid[id2[k]] + mid[id2[k + 1]]) / 2
      reg_ev <- reg[findInterval(pos_ev, mid)]
      rec[[length(rec) + 1]] <- data.table(dataset = f, id = colnames(G)[j], chr = ch,
                                           region = reg_ev, ppos = pos_ev,
                                           gap_win = id2[k + 1] - id2[k])
    }
  }
  R <- rbindlist(rec)
  if (!nrow(R)) next
  out[[f]] <- R
  cat("done", f, "\n")
}
CO <- rbindlist(out)
CO[, family := sub("_WATDE.*$", "", dataset)]
CH <- CO[, .(co = .N / uniqueN(id)), by = .(family, chr)][, co := round(co, 2)]
RG <- CO[, .(co = .N / uniqueN(id)), by = .(family, region)][, co := round(co, 2)]
fwrite(CH, file.path(TAB, paste0("Table18_CO_by_chr_family", SFX, ".csv")))
fwrite(RG, file.path(TAB, paste0("Table18b_CO_by_region_family", SFX, ".csv")))
fwrite(CO[, .(dataset, id, chr, region, ppos, gap_win)], file.path(TAB, paste0("Table18c_CO_events", SFX, ".csv")))
cat("rows", nrow(CO), " families", uniqueN(CO$dataset), "\n")
print(RG[, .(mean = round(mean(co),1)), by = region][order(-mean)])