## 70_candidate_report.R -- candidate gene report for KO planning
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); DOC <- file.path(ROOT, "Docs")
A <- fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))
A[, gene_Mb := round(gene_start / 1e6, 2)]
A[, dist_peak_Mb := round((gene_start - peak_Mb * 1e6) / 1e6, 2)]
R <- A[recomb == TRUE & !is.na(protein)]
R <- unique(R, by = "gene_id")
## curated function notes
note <- function(p) {
  p <- tolower(p)
  if (grepl("helicase q-like|recql", p)) return("RecQ 家族 3'-5' DNA 解旋酶；参与 DNA 双链断裂修复与同源重组，维持基因组稳定性（拟南芥 RECQL1 缺失导致 DNA 损伤修复缺陷）")
  if (grepl("histone demethylase", p)) return("组蛋白去甲基化酶（LSD1/KDM1 同源）；去除 H3K4 甲基化，调控染色质状态与发育时序")
  if (grepl("su(var)|suvr|histone-lysine n-methyltransferase", p)) return("组蛋白甲基转移酶（SUVR4 类）；催化 H3K27 单甲基化，参与转座子沉默与异染色质形成")
  if (grepl("ubiquitin-like modifier-activating", p)) return("泛素样修饰激活酶（UBA 家族）；激活泛素样蛋白，参与蛋白降解/自噬")
  if (grepl("e3 ubiquitin", p)) return("E3 泛素连接酶；底物特异性泛素化，参与蛋白降解、信号转导与发育")
  if (grepl("polyubiquitin", p)) return("多聚泛素；蛋白酶体降解底物标签，维持蛋白稳态")
  if (grepl("atg12", p)) return("ATG12 泛素样蛋白；自噬通路共价连接组分")
  return("—")
}
R[, effect := sapply(protein, note)]
R[, is_helicase := grepl("helicase", protein, ignore.case = TRUE)]
setorder(R, -is_helicase, evalue)
fwrite(R, file.path(TAB, "Table58_candidate_gene_report.csv"))
cat("candidates:", nrow(R), "\n")
for (i in seq_len(nrow(R))) {
  cat(sprintf("%-28s %-4s %7.2f Mb (peak+-%.2f) | %-5s | %-8s | %s | %s | %.1g\n",
      R$gene_id[i], R$chr[i], R$gene_Mb[i], abs(R$dist_peak_Mb[i]), R$direction[i],
      R$family[i], substr(R$protein[i], 1, 40), R$species[i], R$evalue[i]))
}
