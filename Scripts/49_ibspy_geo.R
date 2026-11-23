## 49_ibspy_geo.R -- haplotype diversity and geographic structure at CO-QTL windows
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
PASS <- fread(file.path(TAB, "Table37_WATDE_country.csv"))
coord <- data.table(
  country = c("India","Spain","China","USSR","Iran","Yugoslavia","Portugal","Afghanistan","Australia",
              "Greece","France","Poland","Morocco","Italy","Turkey","Tunisia","Bulgaria","Canary Islands",
              "Crete","Iraq","Hungary","Romania","Algeria","Egypt","Pakistan","Israel","Cyprus","Syria",
              "Lebanon","Jordan","Saudi Arabia","Ethiopia","Kenya","South Africa","Argentina","Chile",
              "Mexico","USA","Canada","Germany","Netherlands","Belgium","UK","Sweden","Denmark",
              "Finland","Norway","Austria","Switzerland","Czechoslovakia","Yugoslavia","Albania"),
  lat = c(22,40,35,55,32,44,39.5,34,25,39,46,52,32,42,39,34,42.7,28,35,33,47,46,28,27,30,31,35,35,
          34,31,24,9,0,-30,-35,-30,23,39,45,51,52,51,54,55,56,60,62,47.5,47,49.5,44,41),
  lon = c(79,-4,105,40,53,20,-8,66,134,22,2,19,-6,12.5,35,9,25.5,-16,25,44,19,25,3,30,70,35,33,38,
          36,36,45,40,38,28,25,-71,-102,-98,-106,10,5,4,-3,15,10,25,10,47,8,17.5,20,20)
)
PASS <- merge(PASS, coord, by = "country", all.x = TRUE)
cat("accessions with coords:", sum(!is.na(PASS$lat)), "of", nrow(PASS), "\n")
## CO-QTL windows: 20 Mb centred on each suggestive QTL peak (physical)
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))
MK <- L$markers[!is.na(g_rank), .(dataset, marker, chr, cM = g_rank, ppos)]
MK[, family := sub("_WATDE.*$", "", dataset)]
PK <- fread(file.path(TAB, "Table30b_significant_peaks.csv"))
PK[, phys := NA_real_]
for (i in seq_len(nrow(PK))) {
  d <- MK[family == PK$family[i] & chr == PK$chr[i]]
  if (nrow(d)) PK$phys[i] <- d[which.min(abs(cM - PK$pos[i])), ppos]
}
QTLW <- unique(PK[!is.na(phys), .(chr, start = phys - 10e6, end = phys + 10e6, type = "CO-QTL")])
QTLW[, start := pmax(0, start)]
## random control windows (20 Mb) on the same chromosomes
set.seed(7)
ctrl <- rbindlist(lapply(unique(QTLW$chr), function(ch) {
  mx <- max(MK[chr == ch]$ppos, na.rm = TRUE)
  st <- seq(0, max(mx - 20e6, 0), by = 20e6)
  data.table(chr = ch, start = st, end = st + 20e6, type = "control")
}))
WINS <- rbind(QTLW, ctrl)
metric <- function(g, lat, lon) {
  ok <- !is.na(g) & !is.na(lat)
  g <- g[ok]; lat <- lat[ok]; lon <- lon[ok]
  if (length(g) < 100) return(NULL)
  n <- length(g); idx <- seq_len(n)
  set.seed(1); np <- 20000
  a <- sample(idx, np, TRUE); b <- sample(idx, np, TRUE)
  keep <- a != b
  share <- g[a[keep]] == g[b[keep]]
  dlat <- abs(lat[a[keep]] - lat[b[keep]])
  small <- mean(share[dlat <= 5]); large <- mean(share[dlat >= 15])
  list(n = n, n_group = length(unique(g)),
       shannon = { p <- table(g) / n; -sum(p * log(p)) },
       share_small = small, share_large = large,
       structure = if (is.finite(large) && large > 0) small / large else NA_real_)
}
res <- list()
fs <- list.files(WD, pattern = "ibs_haplotypes\\.tsv\\.gz$", full.names = TRUE)
fs <- fs[!grepl("\\(1\\)", fs)]
for (f in fs) {
  ch <- sub("^chr([0-9][ABD]).*$", "\\1", basename(f))
  w <- WINS[chr == ch]
  if (!nrow(w)) next
  d <- fread(f, select = c("start", "end", "query", "group"))
  d <- d[query %in% PASS$wATDE]
  d[, group := as.integer(group)]
  for (i in seq_len(nrow(w))) {
    dd <- d[end >= w$start[i] & start <= w$end[i]]
    if (!nrow(dd)) next
    ## one group per accession per window: majority
    g <- dd[, .(group = names(sort(table(group), decreasing = TRUE))[1]), by = query]
    m <- merge(g, PASS[, .(query = wATDE, lat, lon)], by = "query")
    met <- metric(as.integer(m$group), m$lat, m$lon)
    if (is.null(met)) next
    res[[length(res) + 1]] <- data.table(chr = ch, type = w$type[i],
      start = round(w$start[i] / 1e6, 1), end = round(w$end[i] / 1e6, 1),
      n = met$n, n_group = met$n_group, shannon = round(met$shannon, 3),
      share_small = round(met$share_small, 4), share_large = round(met$share_large, 4),
      structure = round(met$structure, 3))
  }
  cat("done", ch, "\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table38_IBSpy_geo_structure.csv"))
print(R[, .(n_win = .N, shannon = round(mean(shannon), 3), structure = round(mean(structure, na.rm = TRUE), 3)), by = type])
cat("\nWilcoxon shannon:", wilcox.test(shannon ~ type, data = R)$p.value, "\n")
cat("Wilcoxon structure:", wilcox.test(structure ~ type, data = R)$p.value, "\n")