# Tables 索引

## 主结果

| 表 | 文件 | 内容 | 脚本 |
|---|---|---|---|
| Table0 | `Table0_cohort.csv` | 家族入组（>500 标记）与保留标记数 | 23_co_final.R |
| Table1a | `Table1a_markers_{AXIOM.1,AXIOM.2,KASP}.csv` | 标记、共识位置、物理位置、各家系遗传距离 | 03_map_tables.R |
| Table1b | `Table1b_map_statistics.csv` | 每条染色体与基因组的图谱统计 | 03_map_tables.R |
| Table1c | `Table1c_map_audit.csv` | 图谱审计（单调性、tie 比例、cM/Mb） | 00_audit_maps.R |
| Table2a | `Table2a_phenotypes_by_line.csv` | 品系表型均值 | 04_phenotypes.R |
| Table2b | `Table2b_phenotype_coverage.csv` | 表型覆盖情况 | 04_phenotypes.R |
| Table2c | `Table2c_QTL_family_consensus_vs_original.csv` | 各家族 QTL（共识 vs 原始图谱） | 05c_qtl_projected.R |
| Table2d | `Table2d_QTL_profiles.csv` | QTL LOD profile（置信区间来源） | 05c_qtl_projected.R |
| Table2e | `Table2e_QTL_consensus_vs_original_summary.csv` | 共识/原始 QTL 汇总 | 07_qtl_figures.R |
| Table2f | `Table2f_metaQTL.csv` | Meta-QTL 聚类 | 09_meta_analysis.R |
| Table2g | `Table2g_QTL_per_chromosome.csv` | 每染色体 QTL 计数 | 09_meta_analysis.R |
| Table3 | `Table3_marker_coverage.csv` | 标记覆盖统计 | 03_map_tables.R |
| Table4 | `Table4_map_length_by_family.csv` | 家族图谱长度 | 10_outlier_family.R |
| Table5 | `Table5_genotype_quality.csv` | 家族基因型质量 | 10_outlier_family.R |
| Table6 | `Table6_chromosome_ends.csv` | 染色体末端统计 | 03_map_tables.R |
| Table7 | `Table7_LD_decay.csv` | LD 衰减（1 Mb 分箱） | 12_ld_decay.R |
| Table8 | `Table8_window_from_LD.csv` | 由 LD 确定的窗口长度（r²=0.8/0.6） | 12_ld_decay.R |
| Table17 | `Table17_CO_direct_vs_window_{15Mb,20Mb}.csv` | 直接法 vs 窗口法（敏感性分析） | 23_co_final.R |
| Table18 | `Table18_CO_by_chr_family.csv` | 家族 × 染色体 CO | 26_co_by_chr_region.R |
| Table18b | `Table18b_CO_by_region_family.csv` | 家族 × 区域 CO | 26_co_by_chr_region.R |
| Table18c | `Table18c_CO_events.csv` | CO 事件（位置、区域、跨度） | 26_co_by_chr_region.R |
| Table19 | `Table19_PD_by_chr.csv` / `Table19b_PD_by_region.csv` | 亲本多态性 PD | 28_pd.R |
| Table20 | `Table20_CO_QTL_profiles.csv` | 四类 CO 表型 QTL profile | 32_co_qtl.R |
| Table20b | `Table20b_CO_QTL_peaks.csv` | 四类 CO 表型 QTL 峰 | 32_co_qtl.R |
| Table21 | `Table21_allele_effects.csv` | 等位频率与等位效应 | 34_allele_effect.R |
| Table22 | `Table22_CO_region_rate.csv` | 区域 CO（速率积分，主结果） | 35_co_rate.R |
| Table22b | `Table22b_CO_10Mb_rate.csv` | 10 Mb 分箱 CO 密度 | 35_co_rate.R |
| Table22c | `Table22c_CO_rate_individual.csv` | 个体级区域 CO（QTL 表型） | 35_co_rate.R |
| Table23 | `Table23_informative_windows.csv` | 每个体可用纯合窗口数（QTL 协变量） | 37_window_covariate.R |

## 归档

`_archive/`：早期被取代的 CO 统计版本（逐标记计数、5.5 Mb 窗口、5 Mb 过渡版）、旧版染色体汇总等，仅作追溯，不用于论文。

## 冻结后新增（2026-11-13）

| 表 | 文件 | 内容 |
|---|---|---|
| Table17 | `Table17_CO_direct_vs_window_15Mb.csv` | 直接法 vs 窗口法（15 Mb，主分析） |
| Table24 | `Table24_PD_tree_bootstrap.csv` | PD 树各节点自展支持率 |
| Table25 | `Table25_PD_clusters.csv` | 家族 → 两大簇归属 |
| Table26 | `Table26_cluster_CO_outlier_test.csv` | 簇间 CO 差异（原始/剔除离群） |
| Table27 | `Table27_cluster_CO_ancova.csv` | ANCOVA（控制窗口数/标记数） |
| Table28 | `Table28_CO_rarefied.csv` | 稀释到 700 标记后的个体 CO |
| Table29 | `Table29_CO_by_family_rarefied.csv` | 稀释后家族均值 |

## 补充分析表（2026-11-13/14）

| 表 | 文件 | 内容 | 脚本 |
|---|---|---|---|
| Table30 | `Table30_CO_metaQTL.csv` | CO 四类表型 meta-QTL 聚类 | 40_co_meta.R |
| Table30b | `Table30b_CO_QTL_suggestive_peaks.csv` | 提示性 CO-QTL 峰（LOD ≥ 2） | 40_co_meta.R |
| Table31 | `Table31_CO_QTL_fitness_allele_effects.csv` | CO-QTL 等位基因对 fitness 的家族内效应 | 42_co_qtl_fitness.R |
| Table32 | `Table32_CO_fitness_correlation.csv` / `Table32b` | 个体 CO–fitness 相关及 Fisher 合并 | 43_co_fitness_correlation.R |
| Table33 | `Table33_IBSpy_window_diversity.csv` | IBSpy 窗口多样性 | 44_ibspy_diversity.R |
| Table34 | `Table34_CO_metaQTL_physical.csv` / `Table34b` | meta-QTL 物理区间与多样性检验 | 45_ibspy_enrichment.R |
| Table35 | `Table35_published_QTL_cohort.csv` | 已发表 QTL（18 家族，1,695 个） | 46_published_qtl.R |
| Table36 | `Table36_CO_QTL_coloc_published.csv` | 共定位二项检验 | 47_co_qtl_published_coloc.R |
| Table37 | `Table37_WATDE_country.csv` | WATDE → 国别 | 49_ibspy_geo.R |
| Table38 | `Table38_IBSpy_geo_structure.csv` | 窗口多样性/地理结构 | 49_ibspy_geo.R |
| Table39 | `Table39_CO_QTL_trait_scan.csv` | CO-QTL 与 24 性状富集扫描 | 50_co_qtl_trait_scan.R |
| Table40 | `Table40_precise_coord.csv` | 精确坐标材料 | 51_ibspy_precise_geo.R |
| Table41 | `Table41_precise_geo_mantel.csv` | 精确坐标 Mantel 检验 | 51_ibspy_precise_geo.R |
| Table42 | `Table42_RIL_phenotypes_by_season.csv` | RIL 分季节表型（抽穗期/株高/产量/千粒重） | 52_ril_heading_date.R |
| Table43 | `Table43_RIL_phenotypes_multiseason.csv` | RIL 多年平均表型 | 52_ril_heading_date.R |

## 候选基因、LD 与传递分析（2026-11-15）

| 表 | 文件 | 内容 |
|---|---|---|
| Table56 | `Table56_significant_QTL_genes_5Mb.csv` | 9 个显著 QTL ±5 Mb 内 620 个基因 |
| Table57 | `Table57_significant_QTL_annotation.csv` | 上述基因的 Swiss-Prot 注释 |
| Table58 | `Table58_candidate_gene_report.csv` | 14 个候选基因（含功能与方向） |
| Table59 | `Table59_RECQL1_hits.tsv` | RECQL1 的 RECQ 家族比对 |
| Table60/60b | `Table60b_homeologs.tsv` | 候选基因在 Paragon 蛋白组中的同源拷贝 |
| Table61 | `Table61_homeolog_gene_models.csv` | 六个拷贝的基因模型 |
| Table62/62b | `Table62*_RECQL1_*.csv` | RECQL1 4A/4D vs 4B 氨基酸差异 |
| Table63/63b | `Table63*_LSD1_*.csv` | LSD1 7A/7B vs 7D 氨基酸差异 |
| Table64 | `Table64_Watkins_QTL_4A_candidates.csv` | Watkins 增效 QTL 的 18 个基因 |
| Table65 | `Table65_LSD1_hits.tsv` | LSD1 的 KDM1 家族比对 |
| Table66 | `Table66_7D_LOD_Hd_covariate.csv` | 7D 加抽穗期协变量的 LOD 曲线 |
| Table67/68 | `Table67/68_*` | 全基因组 LD / 块（按 QTL 基因型） |
| Table69 | `Table69_QTL_allele_transmission.csv` | 9 个 QTL 的等位传递检验 |
| Table70 | `Table70_double_high_blocks.csv` | 高重组等位数 vs 块统计 |
| Table73/74 | `Table73/74_*` | QTL 染色体局部 LD / 块 |
| Table75 | `Table75_cross_family_replication.csv` | 跨家系复现可行性 |
| Table76/77/78 | `Table76/77/78_*` | 家系内 LD / 块 / 检验 |
| Table79/80 | `Table79/80_*` | LD 衰减距离（局部 / 全基因组） |
| Table81/82 | `Table81/82_*` | 残余杂合（总体 / 分区域） |
