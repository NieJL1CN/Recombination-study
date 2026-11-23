# Scripts — 运行顺序与产出对照

所有脚本用 R 4.5.1 运行，工作目录为项目根目录：

```powershell
cd D:\CAAS\Recombination\WatkinsRILsRecombinationStudy
& 'D:\ProgramFiles\R-4.5.1\bin\Rscript.exe' 'Scripts\01_build_maps.R'
```

## 主流程（按编号顺序执行）

| 脚本 | 作用 | 主要产出 |
|---|---|---|
| `00_audit_maps.R` | 共识图谱质量审计（单调性、tie 比例、cM/Mb） | `Tables/Table1c_map_audit.csv` |
| `01_build_maps.R` | 构建共识图谱：物理排序 → 面板内秩匹配 → LPmerge → 物理单调化 → 归零 | `Source/Derived/maps_min1.rds`、`consensus_maps_min1.csv` |
| `02_map_figures.R` | 图 1a（密度）、图 1b（连锁图）、图 1d（共线性 3×7） | `Figs/Fig1a/1b/1d_*.png` |
| `02b_cmplot.R` | CMplot 密度图（物理坐标 + 5 Mb 分箱） | `Figs/Fig1a_cmplot_*.png` |
| `03_map_tables.R` | 表 1a（标记/共识位置/物理位置/各家系遗传距离）、表 1b（统计） | `Tables/Table1a_*.csv`、`Table1b_map_statistics.csv` |
| `04_phenotypes.R` | 公共表型整理 | `Tables/Table2a/2b_*.csv` |
| `05_qtl_family.R` | 家系 QTL（原始图谱） | `Tables/Table2c_*.csv` |
| `05c_qtl_projected.R` | 家系 QTL（投影到共识图谱）与原始图谱对比 | `Tables/Table2c_QTL_family_consensus_vs_original.csv`、`Table2d_QTL_profiles.csv` |
| `06_linkage_heatmap.R` | 图 1c（连锁图胶囊热图） | `Figs/Fig1c_linkagemap_*.png` |
| `07_qtl_figures.R` | 图 2a（曼哈顿，Y 轴统一） | `Figs/Fig2a_manhattan_*.png` |
| `08_mpp_qtl.R` | statgenMPP 多亲本 QTL（样本量受限，附录） | `Tables/Table2h_*.csv` |
| `09_meta_analysis.R` | Meta-QTL 聚类 | `Tables/Table2f_metaQTL.csv`、`Table2g_*.csv` |
| `10_outlier_family.R` | 图谱长度/基因型质量异常家系 | `Tables/Table4/Table5_*.csv`、`Figs/Fig1e_*.png` |
| `11_meta_map.R` | **图 2c**：平行染色体比较图（白条+黑边、标记横线、红色 QTL 区间、同源标记连线、Y 倒置、theme_void） | `Figs/Fig2c_metaQTL_<chr>.png` |
| `12_ld_decay.R` | LD 衰减（1 Mb 分箱） | `Tables/Table7_LD_decay.csv`、`Table8_window_from_LD.csv` |
| `17_ld_curve.R` | 图 3f（LD 曲线，标注 r²=0.8/0.6） | `Figs/Fig3f_LD_decay.png` |
| `23_co_final.R` | CO 直接计数法 vs 窗口法对比 | `Tables/Table17*_CO_direct_vs_window_*.csv` |
| `26_co_by_chr_region.R` | 窗口法统计 CO 事件（按染色体/区域） | `Tables/Table18*_CO_*.csv` |
| `35_co_rate.R` | **核心**：参考方法速率积分（每 CO 贡献 1/跨度，按 Mb 铺开） | `Tables/Table22*_CO_rate*.csv` |
| `36_rate_plot.R` | **图 3d**：全基因组 CO 密度（1A 在下、7D 在上，窄版） | `Figs/Fig3d_CO_genome_distribution.png` |
| `27_co_heatmap.R` | 图 3c：染色体 × 家族 CO 热图（格内标值） | `Figs/Fig3c_CO_heatmap.png` |
| `37_window_covariate.R` | 每个体可用纯合窗口数（QTL 协变量） | `Tables/Table23_informative_windows.csv` |
| `32_co_qtl.R` | **图 3g**：四类 CO 表型 QTL（远端/中部/着丝粒/全基因组，含窗口数协变量） | `Tables/Table20*_CO_QTL_*.csv` |
| `33_co_qtl_figures.R` | 四类 QTL 曼哈顿图（Y 轴统一） | `Figs/Fig3i_CO_QTL_*.png` |
| `34_allele_effect.R` | 等位基因频率、偏分离、等位效应检查 | `Tables/Table21_allele_effects.csv`、`Figs/Fig4a/4b_*.png` |
| `28_pd.R`、`30_pd_tree.R` | 亲本多态性 PD（按染色体/区域）、PD 热图与聚类树 | `Tables/Table19*.csv`、`Figs/Fig3f_PD_*.png` |
| `31_regression.R` | CO ~ 染色体长度、CO ~ PD 回归 | `Figs/Fig3e_CO_vs_*.png` |

## 关键环境变量

| 变量 | 默认 | 说明 |
|---|---|---|
| `CO_WIN` | `5` | CO 窗口长度（Mb） |
| `SUFFIX` | 空 | 输出文件后缀（如 `_5Mb`） |
| `MINFAM` | `1` | LPmerge 每个共识位置的家族数下限 |

## 归档

`_archive/` 存放早期被取代的脚本（5.5 Mb 窗口、逐标记 CO、错位的中间版本），仅作追溯，不参与最终结果。

## 新增脚本（2026-11-13/14，补充分析）

| 脚本 | 作用 | 主要产出 |
|---|---|---|
| `40_co_meta.R` | CO 四类表型的 meta-QTL 聚类（≥2 家族） | `Table30/30b`（只出表，图由 41 绘制） |
| `41_co_meta_map.R` | **Fig3i**：CO QTL + meta-QTL 簇图（区间=峰标记±1 标记） | `Fig3i_CO_metaQTL.png` |
| `42_co_qtl_fitness.R` | CO-QTL 等位基因对 PH/GY/GW 的家族内效应 | `Table31` |
| `43_co_fitness_correlation.R` | 个体 CO 与 fitness 性状的家族内相关 | `Table32/32b`、`Fig4c` |
| `44_ibspy_diversity.R` | IBSpy 窗口级单倍型多样性 | `Table33` |
| `45_ibspy_enrichment.R` | 20 Mb 窗 CO-QTL vs 背景多样性检验 | `Table34/34b`、`Fig4d` |
| `46_published_qtl.R` | 提取 18 家族的已发表 QTL（Cheng et al. 2024 S17） | `Table35` |
| `47_co_qtl_published_coloc.R` | CO-QTL 与已发表 QTL 共定位（单性状版） | `Table36` |
| `48_meta_map_hd.R` | **Fig2c**：性状 QTL 图（GY/GW + 已发表 HD，区间=峰标记±1 标记） | `Fig2c_metaQTL_map.png` |
| `49_ibspy_geo.R` | IBSpy 单倍型多样性/地理结构（国别经纬度） | `Table37/38` |
| `50_co_qtl_trait_scan.R` | CO-QTL 与 24 个性状的全扫描 | `Table39`、`Fig4e` |
| `52_ril_heading_date.R` | 从 JIC 多年份原始表提取 RIL 抽穗期/株高/产量/千粒重 | `Table42/43` |
| `51_ibspy_precise_geo.R` | 精确坐标子集（37 材料）Mantel 检验 | `Table40/41` |

**注**：`11_meta_map.R` 已归档（`_archive/11_meta_map_OLD.R`），Fig2c 由 `48_meta_map_hd.R` 出图。

## 候选基因与 LD 分析脚本（2026-11-15）

| 脚本 | 作用 |
|---|---|
| `54_paragon_candidates.R` / `57_paragon_candidates_tight.R` | 从 Paragon GFF 提取 QTL 区间基因（±5 Mb / ±1 Mb） |
| `58/62_export_candidates.R` | 导出候选蛋白 FASTA |
| `60_recombination_genes.R` | UniProt 获取小麦重组/减数分裂参考基因 |
| `61/63/65_annotate*.R` | DIAMOND 注释与关键词筛选 |
| `66-68_mlh1*.R` | MLH1 位点核验（已否定） |
| `69_significant_QTL_genes_5Mb.R` | 9 个显著 QTL 的基因提取 |
| `72_candidate_report_final.R` | 生成候选基因报告 |
| `74_homeolog_structure.R` / `75_homeolog_diff.R` | 三个拷贝的基因模型与氨基酸差异 |
| `77/78_LSD1*.R` | LSD1 亚家族、物理位置、发育关联 |
| `79_7D_covariate_scan.R` | 7D 加抽穗期协变量的区间作图 |
| `80/84/86_*ld_blocks*.R` | LD 与血缘块（全基因组 / 局部 / 家系内） |
| `81_transmission.R` | QTL 等位传递检验 |
| `82_double_high_blocks.R` | 双高组合块统计 |
| `85_cross_family_replication.R` | 跨家系复现可行性 |
| `88_genomewide_ld_decay.R` | 全基因组 LD 衰减曲线与 d50 |
| `89/90_*heterozygosity*.R` | 残余杂合（总体 / 分区域） |
