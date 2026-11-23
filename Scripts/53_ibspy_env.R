## 53_ibspy_env.R -- improved environment-gradient test for CO-QTL windows
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
## passport with country centroids
PASS <- fread(file.path(TAB, "Table37_WATDE_country.csv"))
coord <- data.table(
  country = c("India","Spain","China","USSR","Iran","Yugoslavia","Portugal","Afghanistan","Australia",
              "Greece","France","Poland","Morocco","Italy","Turkey","Tunisia","Bulgaria","Canary Islands",
              "Crete","Iraq","Hungary","Romania","Algeria","Egypt","Pakistan","Israel","Cyprus","Syria",
              "Lebanon","Jordan","Saudi Arabia","Ethiopia","Kenya","South Africa","Argentina","Chile",
              "Mexico","USA","Canada","Germany","Netherlands","Belgium","UK","Sweden","Denmark","Finland",
              "Norway","Austria","Switzerland","Czechoslovakia","Albania"),
  lat = c(22,40,35,55,32,44,39.5,34,25,39,46,52,32,42,39,34,42.7,28,35,33,47,46,28,27,30,31,35,35,
          34,31,24,9,0,-30,-35,-30,23,39,45,51,52,51,54,55,56,60,62,47.5,47,49.5,41),
  lon = c(79,-4,105,40,53,20,-8,66,134,22,2,19,-6,12.5,35,9,25.5,-16,25,44,19,25,3,30,70,35,33,38,
          36,36,45,40,38,28,25,-71,-102,-98,-106,10,5,4,-3,15,10,25,10,47,8,17.5,20)
)
PASS <- merge(PASS, coord, by = "country", all.x = TRUE)
PASS <- PASS[!is.na(lat)]
## windows
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
## pairwise geographic distance among accessions
gcd <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
ACC <- unique(PASS, by = "wATDE")
lat <- ACC$lat; lon <- ACC$lon
Dg <- outer(seq_along(lat), seq_along(lat), Vectorize(function(i, j) gcd(lat[i], lon[i], lat[j], lon[j])))
res <- list()
fs <- list.files(WD, pattern = "ibs_haplotypes\\.tsv\\.gz$", full.names = TRUE)
fs <- fs[!grepl("\\(1\\)", fs)]
for (f in fs) {
  ch <- sub("^chr([0-9][ABD]).*$", "\\1", basename(f))
  w <- WINS[chr == ch]
  if (!nrow(w)) next
  d <- fread(f, select = c("start","end","query","group"))
  d <- d[query %in% ACC$wATDE]
  d[, group := as.integer(group)]
  for (i in seq_len(nrow(w))) {
    dd <- d[end >= w$start[i] & start <= w$end[i]]
    if (!nrow(dd)) next
    g <- dd[, .(group = names(sort(table(group), decreasing = TRUE))[1]), by = query]
    m <- merge(g, ACC[, .(query = wATDE, lat)], by = "query")
    if (nrow(m) < 200) next
    idx <- match(m$query, ACC$wATDE)
    Dm <- Dg[idx, idx]
    dgm <- as.numeric(Dm[upper.tri(Dm)])
    gv <- as.integer(m$group)
    S <- outer(gv, gv, "==")
    sg <- as.numeric(S[upper.tri(S)])
    mr <- suppressWarnings(cor(sg, dgm))
    ## latitudinal differentiation: Fst-like ratio
    lv <- m$lat
    within <- outer(lv, lv, function(a, b) abs(a - b) <= 5)[upper.tri(S)]
    far <- outer(lv, lv, function(a, b) abs(a - b) >= 15)[upper.tri(S)]
    hs <- mean(sg); hw <- mean(sg[within]); hf <- mean(sg[far])
    Fst <- if (is.finite(hs) && hs > 0) (hs - mean(c(hw, hf))) / hs else NA_real_
    ## latitude slope of the dominant haplotype frequency
    top <- as.integer(names(sort(table(gv), decreasing = TRUE))[1])
    bin <- as.integer(gv == top)
    fit <- suppressWarnings(glm(bin ~ lv, family = binomial()))
    slope <- unname(coef(fit)[2])
    res[[length(res) + 1]] <- data.table(chr = ch, type = w$type[i],
      start = round(w$start[i] / 1e6, 1), end = round(w$end[i] / 1e6, 1),
      n = nrow(m), n_group = uniqueN(gv),
      shannon = { p <- table(gv) / nrow(m); round(-sum(p * log(p)), 3) },
      mantel_r = round(mr, 3), Fst_lat = round(Fst, 4), lat_slope = round(slope, 4))
  }
  cat("done", ch, "\n")
}
R <- rbindlist(res)
fwrite(R, file.path(TAB, "Table44_IBSpy_env_metrics.csv"))
print(R[, .(n_win = .N, n_acc = round(mean(n)),
            shannon = round(mean(shannon), 3),
            mantel = round(mean(mantel_r, na.rm = TRUE), 3),
            Fst_lat = round(mean(Fst_lat, na.rm = TRUE), 4),
            lat_slope = round(mean(abs(lat_slope), na.rm = TRUE), 4)), by = type])
for (v in c("shannon","mantel_r","Fst_lat","lat_slope")) {
  cat(v, "Wilcoxon p =", signif(wilcox.test(get(v) ~ type, data = R)$p.value, 3), "\n")
}