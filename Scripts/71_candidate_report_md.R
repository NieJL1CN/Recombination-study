## 71_candidate_report_md.R -- full candidate report (CSV + Markdown)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); DOC <- file.path(ROOT, "Docs")
A <- fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))
H <- fread(file.path(TAB, "Table56c_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, hit_prot := sub("^[^ ]+ ", "", desc)]
H[, hit_sym := sub(" .*$", "", hit_prot)]
H[, hit_species := sub(".*OS=([^ ]+) OX=.*", "\\1", desc)]
H[, hit_prot := sub(" OS=.*$", "", hit_prot)]
setorder(H, gene_id, -bits)
H1 <- H[, head(.SD, 1), by = gene_id]
R <- merge(A[recomb == TRUE & !is.na(protein)], H1[, .(gene_id, hit_sym, hit_species, hit_prot, pident, alen)], by = "gene_id")
R <- unique(R, by = "gene_id")
R[, gene_Mb := round(gene_start/1e6, 2)]
R[, dist_Mb := round((gene_start - peak_Mb*1e6)/1e6, 2)]
note <- function(p) {
  p <- tolower(p)
  if (grepl("helicase q-like|recql", p)) return("RecQ-family 3'-5' DNA helicase; DSB repair and homologous recombination; loss causes DNA-damage sensitivity")
  if (grepl("histone demethylase", p)) return("Histone demethylase (LSD1/KDM1-like); removes H3K4 methylation, alters chromatin and developmental timing")
  if (grepl("suvr|histone-lysine n-methyltransferase", p)) return("Histone methyltransferase (SUVR4-like); H3K27 monomethylation, transposon silencing, heterochromatin")
  if (grepl("ubiquitin-like modifier-activating", p)) return("Ubiquitin-like modifier activating enzyme (UBA family); protein degradation/autophagy")
  if (grepl("e3 ubiquitin", p)) return("E3 ubiquitin ligase; substrate-specific ubiquitination, protein turnover and signalling")
  if (grepl("polyubiquitin", p)) return("Polyubiquitin; proteasome targeting, protein homeostasis")
  if (grepl("atg12", p)) return("ATG12 ubiquitin-like protein; autophagy conjugation system")
  return("not annotated")
}
R[, effect := sapply(hit_prot, note)]
R[, in_1Mb := abs(dist_Mb) <= 1]
R[, is_hel := grepl("helicase", hit_prot, ignore.case = TRUE)]
setorder(R, -is_hel, evalue)
fwrite(R, file.path(TAB, "Table58_candidate_gene_report.csv"))
con <- file(file.path(DOC, "CANDIDATE_GENES_REPORT.md"), "w", encoding = "UTF-8")
cat("# CO-QTL 候选基因报告（用于基因敲除实验设计）\n\n", file = con)
cat("## 范围与定义\n\n", file = con)
cat("- 仅限 **9 个显著 CO-QTL**（LOD > 500 次置换阈值），窗口为峰位 ±5 Mb\n", file = con)
cat("- 基因坐标来自 Paragon 注释 GCA949126075v1（TraesPARA_EIv1.0_*）\n", file = con)
cat("- 功能注释来自 DIAMOND 比对 Swiss-Prot（同源推断，非实验验证）\n", file = con)
cat("- 方向 = 该 QTL 中增加重组的亲本等位基因（A = Paragon，B = Watkins）\n\n", file = con)
cat("## 候选基因（14 个）\n\n", file = con)
cat("| # | 基因 ID | 染色体 | 位置 (Mb) | 距峰 (Mb) | QTL | 方向 | LOD | 最佳匹配（物种 / 基因） | 一致性 | 功能 |\n", file = con)
cat("|---|---|---|---|---|---|---|---|---|---|---|\n", file = con)
for (i in seq_len(nrow(R))) {
  cat(sprintf("| %d | %s | %s | %.2f | %+.2f | %s-%s | %s | %.2f | %s / %s | %.0f%% | %s |\n",
    i, R$gene_id[i], R$chr[i], R$gene_Mb[i], R$dist_Mb[i], R$trait[i], R$family[i],
    R$direction[i], R$lod[i], R$hit_species[i], R$hit_sym[i], R$pident[i], R$effect[i]), file = con)
}
cat("\n## 敲除优先级建议\n\n", file = con)
cat("1. **TraesPARA_EIv1.0_1336720（RECQL1）** — 4B:168.40 Mb，distal-ParW579-4B（LOD 3.44）。\n", file = con)
cat("   RecQ 解旋酶直接参与同源重组与 DNA 修复，是唯一与重组机制直接相关的候选。建议先单敲 4B 拷贝，再聚合 A/D 拷贝。\n\n", file = con)
cat("2. **TraesPARA_EIv1.0_2066850（SUVR4）** — 6B:462.81 Mb，whole-RobW145-6B（LOD 2.96）。\n", file = con)
cat("   组蛋白 H3K27 甲基化，影响异染色质与交叉分布。\n\n", file = con)
cat("3. **TraesPARA_EIv1.0_2580590（LSD1 类去甲基化酶）** — 7D:393.91 Mb，distal-ParW352-7D（LOD 3.45）。\n", file = con)
cat("   与发育时序相关，符合“重组率与发育时序相关”的假设，但不是直接重组基因。\n\n", file = con)
cat("## 注意事项\n\n", file = con)
cat("- 区间是 ±5 Mb 分析窗口，不是 QTL 统计置信区间（后者仅 ±1 个标记，约 1–2 cM）；\n", file = con)
cat("- 620 个基因中 445 个有命中，仅 14 个命中重组/染色质关键词；未注释基因也可能相关；\n", file = con)
cat("- 敲除前需核对基因模型（可用 IWGSC 同源或 RNA-seq）；\n", file = con)
cat("- 9 个显著 QTL 中 8 个由 Paragon 等位基因增加重组，仅 1 个来自 Watkins；\n", file = con)
cat("- 建议先比较 Paragon 与对应 WATDE 亲本在三拷贝上的序列差异，再决定敲哪个拷贝。\n", file = con)
close(con)
cat("written Docs/CANDIDATE_GENES_REPORT.md\n")