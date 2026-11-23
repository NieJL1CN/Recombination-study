## 76_watkins_4A_report.R -- candidate genes for the Watkins-derived QTL (ParW729-4A)
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); DOC <- file.path(ROOT, "Docs")
A <- fread(file.path(TAB, "Table57_significant_QTL_annotation.csv"))
S <- A[family == "ParW729" & chr == "4A"]
S[, gene_Mb := round(gene_start/1e6, 2)]
S[, dist := round((gene_start - peak_Mb*1e6)/1e6, 2)]
setorder(S, gene_start)
KEEP <- c("gene_id","chr","gene_Mb","dist","family","trait","direction","lod","protein","species","evalue")
fwrite(S[, ..KEEP], file.path(TAB, "Table64_Watkins_QTL_4A_candidates.csv"))
cat("genes:", nrow(S), " with hits:", sum(!is.na(S$protein)),
    " recombination-flagged:", sum(S$recomb), "\n")
cat("peak:", unique(S$peak_Mb), "Mb  LOD:", unique(S$lod), " direction:", unique(S$direction), "\n\n")
for (i in seq_len(nrow(S))) cat(sprintf("%-28s %7.2f Mb (%+.2f) %s\n", S$gene_id[i], S$gene_Mb[i], S$dist[i],
  ifelse(is.na(S$protein), "no Swiss-Prot hit", paste0(S$protein[i], " [", S$species[i], "], E=", signif(S$evalue[i], 2)))))
## append section to the candidate report
con <- file(file.path(DOC, "CANDIDATE_GENES_REPORT.md"), "a", encoding = "UTF-8")
cat("\n## 附：唯一由 Watkins 增效的 QTL（ParW729-4A）\n\n", file = con)
cat("- QTL：distal-ParW729-4A，峰位 530.94 Mb，LOD 3.55，**Watkins 等位增加重组**\n", file = con)
cat("- 窗口内 18 个基因，全部有 Swiss-Prot 命中，**没有一个命中重组/减数分裂关键词**\n", file = con)
cat("- 最值得进一步看的两个：\n", file = con)
cat("  1. **TraesPARA_EIv1.0_1276330（CPL3，RNA Pol II CTD 磷酸酶样 3）** — 转录延伸/染色质调控，可能间接影响 DSB 与 CO 分布；\n", file = con)
cat("  2. **TraesPARA_EIv1.0_1276290（RBR2，视网膜母细胞瘤相关蛋白 2）** — 细胞周期/E2F 调控，减数分裂进入相关；\n", file = con)
cat("- 其余为糖基转移酶、氨肽酶、PIP5K4、PP2C12、ABC 转运蛋白、线粒体外膜转位酶、LURP-one 蛋白及 6 个未注释基因。\n", file = con)
cat("- **结论**：该 QTL 更可能是调控型变异（转录/细胞周期），而非经典重组基因的编码区变异；作为农家种“高重组”等位基因的育种价值仍需功能验证。\n", file = con)
close(con)
cat("\nappended to CANDIDATE_GENES_REPORT.md\n")
