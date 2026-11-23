## 107_extract_vcf_windows.R -- extract variants in candidate gene windows
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]; f <- args[2]
win <- fread(args[3], header = FALSE, col.names = c("gene","start","end"))
samples <- c("paragon","WATDE0075","WATDE0100","Chinese_Spring")
con <- gzfile(f, "r")
idx <- NULL; hdr <- NULL
repeat {
  l <- readLines(con, n = 1, warn = FALSE)
  if (!length(l)) break
  if (startsWith(l, "#CHROM")) { hdr <- strsplit(l, "\t")[[1]]; idx <- match(samples, hdr); break }
}
cat("sample indices:", paste(samples, idx, sep = "=", collapse = ", "), "\n")
out <- list()
repeat {
  ln <- readLines(con, n = 200000, warn = FALSE)
  if (!length(ln)) break
  ok <- startsWith(ln, paste0(chr, "\t"))
  if (!any(ok)) next
  x <- ln[ok]
  pos <- as.integer(sub("^[^\t]+\t([0-9]+)\t.*$", "\\1", x))
  keep <- rep(FALSE, length(pos))
  for (i in seq_len(nrow(win))) keep <- keep | (pos >= win$start[i] & pos <= win$end[i])
  if (!any(keep)) next
  x <- x[keep]; pos <- pos[keep]
  sp <- strsplit(x, "\t", fixed = TRUE)
  rec <- rbindlist(lapply(sp, function(v) {
    g <- v[idx]
    data.table(chr = v[1], pos = as.integer(v[2]), ref = v[4], alt = v[5], info = v[8],
               gt_paragon = sub(":.*$", "", g[1]), gt_w579 = sub(":.*$", "", g[2]),
               gt_w729 = sub(":.*$", "", g[3]), gt_cs = sub(":.*$", "", g[4]))
  }))
  out[[length(out)+1]] <- rec
}
close(con)
R <- rbindlist(out)
R[, gene := ""]
for (i in seq_len(nrow(win))) R[pos >= win$start[i] & pos <= win$end[i], gene := win$gene[i]]
fwrite(R, file.path(TAB, paste0("Table102_VCF_", chr, ".csv")))
cat("variants extracted:", nrow(R), "\n")
if (nrow(R)) print(R[, .N, by = gene])
