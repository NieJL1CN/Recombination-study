## ============================================================================
## 01_build_maps.R -- physical-order marker placement, genetic distances,
## isolation filter (>50 cM both sides), and per-panel LPmerge consensus maps
## Project: WatkinsRILsRecombinationStudy
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(LPmerge)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
GEN  <- file.path(ROOT, "Source/Genotype")
DER  <- file.path(ROOT, "Source/Derived"); dir.create(DER, recursive = TRUE, showWarnings = FALSE)
TAB  <- file.path(ROOT, "Tables"); dir.create(TAB, showWarnings = FALSE)
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))

## ---- 1. marker physical positions from the 35K probe set -------------------
p <- as.data.table(readxl::read_excel(file.path(GEN, "35k_probe_set_IWGSCv1.xlsx"), sheet = 1))
setnames(p, c("snp","chr","pos","affx","seq","type","affy820","ba","iselect","bs_cdb","bs_blast","kasp"))
p[, chr := sub("^chr", "", chr)]
p[, pos := suppressWarnings(as.numeric(pos))]
p <- p[chr %in% CHR & !is.na(pos)]
cand <- rbind(p[, .(id = snp, chr, pos)],
              p[!is.na(bs_cdb) & bs_cdb != "", .(id = bs_cdb, chr, pos)],
              p[!is.na(bs_blast) & bs_blast != "", .(id = bs_blast, chr, pos)],
              p[!is.na(iselect) & iselect != "", .(id = sub("_51$", "", iselect), chr, pos)])
cand[, id := trimws(id)]
LU <- unique(cand[!is.na(id) & id != ""], by = "id")
saveRDS(LU, file.path(DER, "probe_lookup.rds"))
cat("probe lookup ids:", nrow(LU), "\n")

norm_id <- function(x) {
  x <- sub("\\.(1[ABD]|2[ABD]|3[ABD]|4[ABD]|5[ABD]|6[ABD]|7[ABD])$", "", x)
  sub("-(1[ABD]|2[ABD]|3[ABD]|4[ABD]|5[ABD]|6[ABD]|7[ABD])$", "", x)
}

## ---- 2. read all genotype files (map part only) ----------------------------
fs <- sort(list.files(file.path(GEN, "Watkins_RIL_Genotypes"), pattern = "_rqtl_genotypes\\.csv$", full.names = TRUE))
M <- rbindlist(lapply(fs, function(f) {
  con <- file(f, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, n = 3), ",", fixed = TRUE)
  data.table(dataset = sub("_rqtl_genotypes\\.csv$", "", basename(f)),
             marker = hdr[[1]][-1], chrlab = hdr[[2]][-1],
             pos = suppressWarnings(as.numeric(hdr[[3]][-1])))
}))
M[, chr := sub("^([1-7][ABD]).*$", "\\1", chrlab)]
M[, bin := suppressWarnings(as.numeric(sub("^[1-7][ABD]", "", chrlab)))]
M[is.na(bin), bin := 0]
M <- M[!is.na(chr) & chr %in% CHR]
M[, clean := norm_id(marker)]
M[LU, on = .(clean = id), pchr := i.chr]
M[LU, on = .(clean = id), ppos := i.pos]
M <- M[!is.na(ppos) & pchr == chr]
cat("marker rows with physical position:", nrow(M), " datasets:", uniqueN(M$dataset), "\n")

## panel = independent SNP array
M[, panel := fifelse(grepl("AX-", marker) & grepl("\\.[1-7][ABD]$", marker), "AXIOM.1",
            fifelse(grepl("AX-", marker) & grepl("-[1-7][ABD]$", marker), "AXIOM.2", "KASP"))]
print(M[, .(rows = .N, markers = uniqueN(marker), datasets = uniqueN(dataset)), by = panel])

## ---- 3. per-family genetic distance along the physical order ---------------
## each family: markers sorted by physical position, genetic distances pulled
## from its own map; an isolated marker (> 50 cM to both neighbours) is dropped
setorder(M, dataset, ppos)
M[, gpos := { v <- pos
  ## flip reversed maps: if the family map runs against the physical order
  if (.N >= 5 && cor(pos, ppos, method = "spearman") < 0) v <- max(pos) - pos
  v }, by = dataset]
iso_fix <- function(x, thr = 50, max_iter = 5) {
  for (it in seq_len(max_iter)) {
    idx <- which(!is.na(x)); if (length(idx) < 2) break
    pv <- x[idx]; g <- diff(pv)
    bad <- integer(0)
    ## interior: both neighbours further than thr
    if (length(g) >= 2) {
      l <- c(NA, g); r <- c(g, NA)
      bad <- which(!is.na(l) & !is.na(r) & l > thr & r > thr)
    }
    ## chromosome start / end: single neighbour further than thr
    if (length(g) >= 1 && g[1] > thr) bad <- c(bad, 1L)
    if (length(g) >= 1 && g[length(g)] > thr) bad <- c(bad, length(idx))
    bad <- unique(bad)
    if (!length(bad)) break
    x[idx[bad]] <- NA
  }
  x
}
M[, gpos := iso_fix(gpos), by = dataset]
M[, dropped := is.na(gpos) & !is.na(pos)]
cat("isolated markers set to missing:", sum(M$dropped, na.rm = TRUE), "\n")
M[, g_rank := { v <- gpos; ok <- !is.na(v); if (any(ok)) v[ok] <- sort(v[ok]); v }, by = .(dataset, chr)]
M[, g_rank := g_rank - min(g_rank, na.rm = TRUE), by = .(dataset, chr)]

fwrite(M[, .(dataset, panel, chr, marker, ppos, original_pos = pos, gpos, g_rank, dropped)],
       file.path(DER, "family_marker_positions.csv"))

## ---- 4. per-panel consensus maps ------------------------------------------
MINFAM <- as.integer(Sys.getenv("MINFAM", "2"))
cat("min families per marker:", MINFAM, "\n")
## targeted exclusion: ParW049 7A is 1,139 cM end to end (5x the other families)
EXCL <- data.table(dataset = c("ParW049_WATDE0175_watseq_axiom", "ParW007_WATDE0002_watseq_axiom"),
                   chr = c("7A", "3D"))
run_panel_chr <- function(P, ch) {
  d <- M[panel == P & chr == ch & !is.na(g_rank)]
  if (nrow(EXCL)) d <- d[!EXCL, on = .(dataset, chr)]
  if (uniqueN(d$dataset) < 2 || nrow(d) < 20) return(NULL)
  L <- d[, .(len = max(g_rank) - min(g_rank), n = .N), by = dataset]
  med <- median(L$len)
  L[, ok := len >= 0.15 * med & len <= 7 * med & n >= 8]
  d <- merge(d, L[, .(dataset, ok)], by = "dataset")[ok == TRUE]
  if (uniqueN(d$dataset) < 2) return(NULL)
  w <- dcast(d[, .(marker, dataset, g_rank)], marker ~ dataset, value.var = "g_rank")
  cnt <- rowSums(!is.na(as.matrix(w[, -1])))
  wf <- w[cnt >= MINFAM]
  if (nrow(wf) < 10) return(NULL)
  maps <- lapply(names(wf)[-1], function(id) {
    z <- data.table(marker = wf$marker, position = wf[[id]])[!is.na(position)]
    z <- z[!duplicated(marker)]; as.data.frame(z)
  })
  names(maps) <- names(wf)[-1]
  maps <- maps[sapply(maps, nrow) >= 3]
  if (length(maps) < 2) return(NULL)
  ok <- try(capture.output(res <- LPmerge(Maps = maps, max.interval = 1:3)), silent = TRUE)
  if (inherits(ok, "try-error")) return(NULL)
  cons <- data.table(marker = as.character(res[[length(res)]][[1]]),
                     cM = as.numeric(res[[length(res)]][[2]]))
  cons <- cons[order(cM)]
  ## trim leading / trailing fragments that are separated by a >25 cM gap
  ## and account for less than 15% of the chromosome length
  repeat {
    if (nrow(cons) < 5) break
    L <- max(cons$cM) - min(cons$cM); if (L <= 0) break
    g <- diff(cons$cM)
    thr <- 10
    big <- which(g > thr)
    if (!length(big)) break
    i <- big[1]
    if (i <= 5 && cons$cM[i] - min(cons$cM) < 0.10 * L) { cons <- cons[-(1:i)]; next }
    j <- big[length(big)]
    if ((nrow(cons) - j) <= 5 && max(cons$cM) - cons$cM[j + 1] < 0.10 * L) { cons <- cons[1:j]; next }
    break
  }
  ## enforce monotonicity by re-assigning the sorted cM values to the
  ## physical order (rank matching; keeps the cM distribution, no plateaus)
  pp <- unique(M[chr == ch, .(marker, ppos)])
  cons <- merge(cons, pp, by = "marker")
  setorder(cons, ppos)
  cons[, cM := sort(cM)]
  cons[, ppos := NULL]
  cons[, cM := cM - min(cM)]
  cons <- cons[order(cM)]
  cons[, chr := ch]
  cons[, panel := P]
  cons[, n_families := length(maps)]
  cons
}
out <- list()
for (P in unique(M$panel)) for (ch in CHR) {
  x <- run_panel_chr(P, ch)
  if (is.null(x)) { cat("skip", P, ch, "\n"); next }
  out[[paste(P, ch)]] <- x
  cat(sprintf("%-8s %-3s markers=%4d len=%6.1f families=%d\n", P, ch, nrow(x), max(x$cM), x$n_families[1]))
}
CM <- rbindlist(out)
fwrite(CM, file.path(DER, paste0("consensus_maps_min", MINFAM, ".csv")))
saveRDS(list(markers = M, consensus = CM), file.path(DER, paste0("maps_min", MINFAM, ".rds")))
cat("consensus rows:", nrow(CM), " panels:", uniqueN(CM$panel), "\n")