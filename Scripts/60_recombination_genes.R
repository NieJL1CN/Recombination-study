## 60_recombination_genes.R -- wheat recombination/meiosis gene reference
suppressPackageStartupMessages({library(data.table); library(httr); library(jsonlite)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables")
queries <- c(
  'organism_id:4565 AND keyword:"Meiosis"',
  'organism_id:4565 AND keyword:"DNA recombination"',
  'organism_id:4565 AND (protein_name:"RAD51" OR protein_name:"DMC1" OR protein_name:"SPO11")',
  'organism_id:4565 AND (protein_name:"MSH4" OR protein_name:"MSH5" OR protein_name:"MLH1" OR protein_name:"MLH3")',
  'organism_id:4565 AND (protein_name:"HEI10" OR protein_name:"ZYP1" OR protein_name:"ASY1" OR protein_name:"REC8" OR protein_name:"FANCM")',
  'organism_id:4565 AND (protein_name:"TOP2" OR protein_name:"CENH3" OR protein_name:"SMC2" OR protein_name:"MRE11" OR protein_name:"RAD50")'
)
get_tsv <- function(q) {
  u <- paste0("https://rest.uniprot.org/uniprotkb/search?query=", URLencode(q, reserved = TRUE),
              "&fields=accession,gene_names,protein_name,xref_ensemblplants&format=tsv&size=100")
  r <- try(GET(u, timeout(60)), silent = TRUE)
  if (inherits(r, "try-error") || status_code(r) != 200) return(NULL)
  txt <- content(r, "text", encoding = "UTF-8")
  if (nchar(txt) < 50) return(NULL)
  d <- fread(text = txt, sep = "\t", quote = "")
  setnames(d, c("Entry", "gene_names", "protein_name", "ep"), skip_absent = TRUE)
  d
}
res <- list()
for (q in queries) {
  d <- get_tsv(q)
  if (is.null(d) || !nrow(d) || !"ep" %in% names(d)) next
  res[[length(res) + 1]] <- d
  cat("query ok ->", nrow(d), "entries\n")
}
U <- unique(rbindlist(res, fill = TRUE), by = "Entry")
cat("unique UniProt entries:", nrow(U), "\n")
ep <- U[!is.na(ep) & ep != "", .(Entry, gene = sub("\\..*$", "", sub(";.*$", "", ep)))]
cat("gene mappings:", nrow(ep), "\n")
getpos <- function(g) {
  id <- sub("\\..*$", "", g)
  u <- paste0("https://rest.ensembl.org/lookup/id/", id,
              "?content-type=application/json;species=triticum_aestivum")
  r <- try(GET(u, timeout(30)), silent = TRUE)
  if (inherits(r, "try-error") || status_code(r) != 200) return(NULL)
  j <- fromJSON(content(r, "text", encoding = "UTF-8"))
  if (is.null(j$seq_region_name)) return(NULL)
  data.table(gene = id, chr = j$seq_region_name, start = j$start, end = j$end)
}
pos <- rbindlist(lapply(head(unique(ep$gene), 150), getpos), fill = TRUE)
cat("positions retrieved:", nrow(pos), "\n")
REF <- merge(ep, U[, .(Entry, gene_names, protein_name)], by = "Entry")
REF <- merge(REF, pos, by = "gene", all.x = TRUE)
REF <- REF[!is.na(chr)]
fwrite(REF, file.path(TAB, "Table49_wheat_recombination_genes.csv"))
cat("reference table:", nrow(REF), "genes with positions\n")
print(REF[, .N, by = chr][order(-N)][1:8])