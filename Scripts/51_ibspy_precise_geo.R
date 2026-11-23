## 51_ibspy_precise_geo.R -- precise-coordinate subset (S23) haplotype vs distance
suppressPackageStartupMessages({library(readxl); library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
## S23 passport with precise coordinates
s23 <- as.data.table(suppressWarnings(read_excel(file.path(WD, "41586_2024_7682_MOESM6_ESM.xlsx"),
                                                 sheet = "S23", skip = 1)))
setnames(s23, c("sample","i19","i16","storecode","watkins_acc","ccode","country","population",
                "locality","coords","elev","growth","ag"))
s23 <- s23[grepl("^WATDE", sample)]
parse_coord <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "0" || !grepl("[0-9]", x)) return(c(NA, NA))
  m <- regmatches(x, regexec("([0-9.]+)\\s*([NS])[, ]*\\s*([0-9.]+)\\s*([EW])", x))[[1]]
  if (length(m) < 5) return(c(NA, NA))
  lat <- as.numeric(m[2]) * ifelse(m[3] == "S", -1, 1)
  lon <- as.numeric(m[4]) * ifelse(m[5] == "W", -1, 1)
  c(lat, lon)
}
ll <- t(sapply(s23$coords, parse_coord))
s23[, lat := ll[, 1]]; s23[, lon := ll[, 2]]
s23 <- s23[!is.na(lat)]
cat("precise-coordinate accessions:", nrow(s23), "\n")
fwrite(s23[, .(sample, country, locality, lat, lon, ag)], file.path(TAB, "Table40_precise_coord.csv"))
## CO-QTL 20 Mb windows
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
PK[, phys := NA_real_]
for (i in seq_len(nrow(PK))) {
  d <- MK[family == PK$family[i] & chr == PK$chr[i]]
  if (nrow(d)) PK$phys[i] <- d[which.min(abs(cM - PK$pos[i])), ppos]
}
QTLW <- unique(PK[!is.na(phys), .(chr, start = pmax(0, phys - 10e6), end = phys + 10e6, type = "CO-QTL")])
set.seed(7)
CTRL <- rbindlist(lapply(unique(QTLW$chr), function(ch) {
  mx <- max(MK[chr == ch]$ppos, na.rm = TRUE)
  st <- seq(0, max(mx - 20e6, 0), by = 20e6)
  data.table(chr = ch, start = st, end = st + 20e6, type = "control")
}))
WINS <- rbind(QTLW, CTRL)
## great-circle helper
gcd <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
mantel_r <- function(g, lat, lon) {
  ok <- !is.na(g)
  g <- g[ok]; lat <- lat[ok]; lon <- lon[ok]
  n <- length(g)
  if (n < 12) return(NA_real_)
  idx <- expand.grid(i = seq_len(n), j = seq_len(n))
  idx <- idx[idx$i < idx$j, ]
  share <- as.numeric(g[idx$i] == g[idx$j])
  dist <- gcd(lat[idx$i], lon[idx$i], lat[idx$j], lon[idx$j])
  if (sd(share) == 0 || sd(dist) == 0) return(NA_real_)
  cor(share, dist)
}
res <- list()
fs <- list.files(WD, pattern = "ibs_haplotypes\\.tsv\\.gz$", full.names = TRUE)
fs <- fs[!grepl("\\(1\\)", fs)]
for (f in fs) {
  ch <- sub("^chr([0-9][ABD]).*$", "\\1", basename(f))
  w <- WINS[chr == ch]
  if (!nrow(w)) next
  d <- fread(f, select = c("start","end","query","group"))
  d <- d[query %in% s23$sample]
  if (!nrow(d)) next
  for (i in seq_len(nrow(w))) {
    dd <- d[end >= w$start[i] & start <= w$end[i]]
    if (!nrow(dd)) next
    g <- dd[, .(group = names(sort(table(group), decreasing = TRUE))[1]), by = query]
    m <- merge(g, s23[, .(query = sample, lat, lon)], by = "query")
    r <- mantel_r(as.integer(m$group), m$lat, m$lon)
    if (is.na(r)) next
    res[[length(res) + 1]] <- data.table(chr = ch, type = w$type[i],
      start = round(w$start[i] / 1e6, 1), end = round(w$end[i] / 1e6, 1),
      n_acc = nrow(m), n_group = uniqueN(m$group), mantel_r = round(r, 3))
  }
  cat("done", ch, "\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table41_precise_geo_mantel.csv"))
print(R[, .(n_win = .N, n_acc = round(mean(n_acc)), mantel = round(mean(mantel_r, na.rm = TRUE), 3)), by = type])
if (nrow(R[type == "CO-QTL"]) >= 3 && nrow(R[type == "control"]) >= 3)
  cat("\nWilcoxon mantel:", wilcox.test(mantel_r ~ type, data = R)$p.value, "\n")
