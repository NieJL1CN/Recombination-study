suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); DOC <- file.path(ROOT, "Docs")
A <- fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))
H <- fread(file.path(TAB, "Table56c_diamond.tsv"), header = FALSE,
           col.names = c("gene_id","target","desc","pident","alen","evalue","bits"))
H[, entry := tstrsplit(target, "|", fixed = TRUE)[[3]]]
H[, hit_sym := sub("_.*$", "", entry)]
H[, code := sub(".*_", "", entry)]
spmap <- c(ARATH = "Arabidopsis thaliana", ORYSJ = "rice (japonica)", ORYSI = "rice (indica)",
           HORVU = "barley", XENLA = "Xenopus", HUMAN = "human", BACSU = "B. subtilis",
           DANRE = "zebrafish", PSEMX = "Pseudomonas", YEAST = "yeast")
H[, hit_species := ifelse(code %in% names(spmap), spmap[code], code)]
H[, hit_prot := sub(" OS=.*$", "", sub("^[^|]+\\|[^|]+\\|[^ ]+ ", "", desc))]
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
  if (grepl("suvr|histone-lysine n-methyltransferase", p)) return("Histone methyltransferase (SUVR4-like); H3K9me1->H3K9me2/3, heterochromatin and transposon silencing")
  if (grepl("ubiquitin-like modifier-activating", p)) return("Ubiquitin-like modifier activating enzyme (UBA family); protein degradation/autophagy")
  if (grepl("e3 ubiquitin", p)) return("E3 ubiquitin ligase; substrate-specific ubiquitination, protein turnover and signalling")
  if (grepl("polyubiquitin", p)) return("Polyubiquitin; proteasome targeting, protein homeostasis")
  if (grepl("atg12", p)) return("ATG12 ubiquitin-like protein; autophagy conjugation system")
  return("not annotated")
}
R[, effect := sapply(hit_prot, note)]
R[, is_hel := grepl("helicase|recql", hit_prot, ignore.case = TRUE)]
## Robigus (RobW families) has only a scaffold-level assembly and no annotation
R[, validation := ifelse(grepl("^RobW", family),
    "deferred (Robigus = GCA_902810685.1, scaffold-level, no annotation)",
    "feasible (Paragon annotation + WatSeq/VCF)")]
setorder(R, -is_hel, evalue)
fwrite(R, file.path(TAB, "Table58_candidate_gene_report.csv"))
con <- file(file.path(DOC, "CANDIDATE_GENES_REPORT.md"), "w", encoding = "UTF-8")
cat("# CO-QTL 候选基因报告（用于基因敲除实验设计）\n\n", file = con)
cat("## 范围与定义\n\n", file = con)
cat("- 仅限 **9 个显著 CO-QTL**（LOD > 500 次置换阈值），窗口为峰位 ±5 Mb\n", file = con)
cat("- 基因坐标来自 Paragon 注释 GCA949126075v1（TraesPARA_EIv1.0_*）\n", file = con)
cat("- 功能注释来自 DIAMOND 比对 Swiss-Prot（同源推断，非实验验证）\n", file = con)
cat("- 方向 = 该 QTL 中增加重组的亲本等位基因：ParW 家系 A = Paragon，RobW 家系 A = **Robigus**，B = Watkins\n\n", file = con)
cat("## 候选基因（14 个）\n\n", file = con)
cat("| # | 基因 ID | 染色体:位置 (Mb) | 距峰 (Mb) | QTL（表型-家系） | 方向 | LOD | 最佳匹配 | 一致性 | 功能 |\n", file = con)
cat("|---|---|---|---|---|---|---|---|---|---|\n", file = con)
for (i in seq_len(nrow(R))) {
  cat(sprintf("| %d | %s | %s:%.2f | %+.2f | %s-%s | %s | %.2f | %s / %s（%s） | %.0f%% | %s |\n",
    i, R$gene_id[i], R$chr[i], R$gene_Mb[i], R$dist_Mb[i], R$trait[i], R$family[i],
    R$direction[i], R$lod[i], R$hit_species[i], R$hit_sym[i], R$hit_prot[i], R$pident[i], R$effect[i]), file = con)
}
cat("\n## 亲本间验证可行性\n\n", file = con)
cat("- ParW 家系：Paragon 有染色体级注释（GCA949126075v1），Watkins 亲本可用 WatSeq/VCF 比对，可做候选基因等位差异验证。\n", file = con)
cat("- RobW 家系：Robigus 仅有 scaffold 级组装（GCA_902810685.1）且无注释，暂不做亲本间候选基因验证；相关 QTL（RobW145-6B）保留在结果中但不进入敲除优先级。\n\n", file = con)
cat("## 敲除优先级建议\n\n", file = con)
cat("1. **TraesPARA_EIv1.0_1336720（RECQL1）** — 4B:168.40 Mb，distal-ParW579-4B（LOD 3.44，Paragon 增效）。\n", file = con)
cat("   RecQ 解旋酶在减数分裂中多为**抗交换因子**（解开 D-loop，抑制非干涉型交换）；Paragon 等位基因增加 CO，提示 Paragon 拷贝可能是功能减弱（LoF）型。\n", file = con)
cat("   注释核验：该蛋白对 Swiss-Prot 的最佳命中为拟南芥 RECQL1（bitscore 684，60% 一致性），其次 RECQL4A（375）、BLM 同源（372–375）、RECQL3（338）、RECQL2（331），支持其为 RECQL1 而非 RECQ4。\n\n", file = con)
cat("2. ~~TraesPARA_EIv1.0_2066850（SUVR4）~~ — 6B:462.81 Mb，whole-RobW145-6B（Robigus 增效）。**暂缓**：Robigus 无注释，无法做亲本间序列验证。\n", file = con)
cat("   植物 SUVR4 主要为 H3K9 甲基转移酶（H3K9me1→me2/3），维持异染色质与转座子沉默；其功能减弱会松动近着丝粒重组抑制，导致全染色体尺度 CO 上升——与该 QTL 的 whole 表型吻合。\n\n", file = con)
cat("3. **TraesPARA_EIv1.0_2580590（LSD1 类）** — 7D:393.91 Mb，distal-ParW352-7D（LOD 3.45）。与发育时序相关，符合“重组率与发育时序相关”的假设，但非直接重组基因。\n\n", file = con)
cat("## 注意事项\n\n", file = con)
cat("- 区间为 ±5 Mb 分析窗口，不是统计置信区间（后者约 1–2 cM）；\n", file = con)
cat("- 620 个基因中 445 个有 Swiss-Prot 命中，仅 14 个命中重组/染色质关键词；未注释基因也可能相关；\n", file = con)
cat("- 敲除前需核对基因模型（IWGSC 同源或 RNA-seq）；\n", file = con)
cat("- 9 个显著 QTL：7 个由 Paragon 增效、1 个由 Robigus 增效（RobW145-6B）、1 个由 Watkins 增效（ParW729-4A）；\n", file = con)
cat("- 建议先比较 Paragon 与对应 WATDE 亲本在三拷贝上的序列差异。\n", file = con)
close(con)
cat("done\n")