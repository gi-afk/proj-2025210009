# 数据卡：Yeast ARS 500 bp windows（v1.0）

本文件说明仓库中数据集的**来源、版本、许可、处理流程与已知边界**。原始文件不随仓库分发，由 `scripts/download_data.py` 下载，并在 `data/raw/download_manifest.json` 记录 URL、字节数、SHA-256 与下载日期。

---

## 1. 数据集概览

| 项目 | 内容 |
|---|---|
| 数据集名称 | Yeast ARS 500 bp window dataset |
| 版本 | v1.0（生成于 2026-11-07） |
| 任务 | 二分类：实验验证的 ARS 窗口 vs 背景窗口 |
| 正样本 | 408 个 500 bp 窗口（OriDB Confirmed ARS 中点） |
| 负样本 A | 408 个全基因组随机窗口 |
| 负样本 B | 408 个距离匹配窗口（与正样本到最近着丝粒/端粒的距离匹配） |
| 染色体 | 1–16（核基因组）；线粒体不参与 |
| 窗口大小 | 500 bp |
| 总窗口数 | 1,224（1:1:1） |
| 核基因组长度 | 12,071,326 bp |
| 随机种子 | 2024 |
| 坐标体系 | R64 / sacCer3（S288C），0-based half-open（内部） |
| 仓库内示例 | `data/sample/sample_dataset.csv`（15 + 15 + 15，仅用于冒烟测试） |

**没有编造、平移或合成任何坐标。** 所有正样本坐标直接来自 OriDB 公开表；所有负样本坐标由脚本在公开参考基因组上按固定规则采样。

---

## 2. 数据来源

### 2.1 参考基因组

| 字段 | 值 |
|---|---|
| 来源 | NCBI RefSeq |
| 组装 | **GCF_000146045.2**（Assembly name: R64；synonym: sacCer3） |
| 物种/菌株 | *Saccharomyces cerevisiae* S288C（Taxid 559292） |
| 组装报告日期 | 2014-12-17 |
| 文件 | `GCF_000146045.2_R64_genomic.fna.gz` |
| URL | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz` |
| 用途 | 提取 500 bp 序列 |
| 许可 | 美国联邦政府资助的公开数据；使用请引用 RefSeq/SGD |

> 说明：最初计划使用 Ensembl Fungi release-59 的 `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz`，但当前环境的 TLS 证书校验失败。NCBI RefSeq 的组装报告明确标注 synonym 为 **sacCer3**，与 Ensembl 的 R64-1-1 使用同一 S288C 参考坐标；因此改用 NCBI。坐标一致性通过正样本 ACS 富集检验（见 §5）间接验证。

染色体命名映射：NCBI 头部为 `chromosome I` … `chromosome XVI`，脚本映射为 `"1"` … `"16"`；线粒体命名为 `MT`，不参与建模。

### 2.2 ARS 坐标

| 字段 | 值 |
|---|---|
| 来源 | OriDB（S. cerevisiae） |
| 端点 | `http://cerevisiae.oridb.org/data_output.php` |
| 参数 | `main=sc_ori&table=sc_ori&format=tab` |
| 访问日期 | 2026-11-07 |
| 站点版本 | 端点未提供版本字符串；以访问日期标识 |
| 总记录 | 829 条（Confirmed 410、Likely 216、Dubious 203） |
| 使用子集 | **status == Confirmed**，410 条 |
| 文件 | `oridb_sc_ori.tsv`（从 HTML `<textarea>` 中提取的 TSV） |
| 许可 | OriDB 面向学术研究免费提供；再分发/商用前请核对官网条款 |

字段：`chr, start, end, name, othernames, status`。坐标按 **1-based、闭区间**发布；脚本转换为中心点 `mid = (start + end) // 2`，窗口为 `[mid - 250, mid + 250)`（0-based half-open）。

交叉核对：`SGD_features.tab` 中有 **352** 条 `ARS` 与 **196** 条 `ARS_consensus_sequence`（覆盖 16 条染色体），与 OriDB 的 Confirmed 集合部分重叠但并非同一份 curated 列表；本项目**不使用 SGD 的 ARS 作为正样本**，只用其着丝粒坐标。两者数量差异是已知的策展差异，不是数据错误。

### 2.3 着丝粒与端粒

| 字段 | 值 |
|---|---|
| 来源 | Saccharomyces Genome Database (SGD) |
| 文件 | `SGD_features.tab` |
| URL | `https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab` |
| 访问日期 | 2026-11-07 |
| 用途 | 16 个着丝粒（`feature type == centromere`），中点 `(start+stop)//2` |
| 端粒 | 以染色体两端（位置 0 与染色体长度）近似；未使用 SGD 的 `telomere` 注释 |
| 许可 | SGD 数据公开；使用请引用 SGD 并遵循其条款 |

端粒近似是本数据集的明确简化：真实端粒重复阵列长度与注释边界会带来数百 bp 到数 kb 的不确定度。此误差只影响距离特征，且对正负样本以同一规则计算。

---

## 3. 处理流程

1. `scripts/download_data.py` 下载基因组、SGD 特征表与 OriDB 表；写 manifest。
2. `scripts/build_dataset.py`：
   - 解析 FASTA，得到 16 条核染色体序列；
   - 读取 OriDB Confirmed 记录，计算中点；
   - 合并中点相距 < 500 bp 的条目（窗口重叠），得到 **408** 个正窗口；
   - 生成 408 个随机负窗口：按染色体长度加权均匀采样，排除距任何正样本中点 < 1,000 bp 的窗口、N 占比 > 2% 的窗口；
   - 生成 408 个距离匹配负窗口：每个正样本在同染色体上寻找一个候选窗口，使其到最近着丝粒/端粒的**实际距离**落在目标距离的容差内（2% → 5% → 10% → 25% → 50%，逐步放宽）；同样排除正样本 1 kb 邻域与 N 多的窗口；
   - 计算 5 个特征：GC、AT、ACS PWM 最大得分、3-mer 熵、log1p 距离；
   - 写 `data/processed/dataset.csv`、`dataset_summary.json` 与 45 窗口示例。
3. `scripts/run_all.py`：按染色体 LeaveOneGroupOut 训练/评估，bootstrap，出图。

---

## 4. 标签与负样本定义

| 类别 | 定义 | 数量 |
|---|---|---|
| `positive` | OriDB `Confirmed` ARS 注释中点 ± 250 bp | 408 |
| `random` | 全基因组均匀随机 500 bp 窗口，距正样本 ≥ 1 kb | 408 |
| `distance_matched` | 每个正样本一个同染色体窗口，最近地标距离匹配 | 408 |

**为什么有两组负样本**：随机负样本可能位于与 ARS 完全不同的基因组环境（GC、基因密度、位置），使模型轻易获得高 AUC；距离匹配负样本固定位置因素，用来检验模型是否真的依赖序列。核心对比见 README §4.3。

---

## 5. 质量检查

| 检查 | 结果 |
|---|---|
| 正样本 ACS 富集（精确 11 bp 共识） | 含 ≥1 个 ACS 的正窗口 **31%**，随机窗口 **2%**（约 15 倍富集） |
| 坐标一致性 | 上述富集只有在坐标与 R64 一致时才成立；结果支持坐标匹配 |
| 正样本 GC 中位数 | 0.328；随机 0.384；匹配 0.386 |
| 正样本 ACS PWM 中位数 | 12.25；随机/匹配 5.63 |
| N 碱基 | 所有入选窗口 N 占比 ≤ 2%；实际最大值为 0 |
| 距离匹配 | 中位数 88,328 vs 87,775 bp；|Δlog1p| p95 = 0.020；<10 kb 比例两边 10.5%；MWU p = 0.997 |
| 染色体验证 | 16 条染色体每条都有正、负样本；交叉验证按染色体留出 |

---

## 6. 已知偏差与限制

- **确认偏差（ascertainment bias）**：只有被实验验证过的 ARS 才会进入正样本；低效、晚期 firing 或重复序列附近的起始点可能系统性缺失。
- **标签宽度**：OriDB 注释区间宽窄不一；统一取中点 ± 250 bp，可能截断较宽的 ARS 或包含邻近序列。
- **负样本不是实验阴性**：随机与距离匹配窗口都只是"背景"，不是"经实验证明不启动复制的片段"。
- **GC/AT 共线**：特征表中二者互补，模型解释需合并。
- **端粒近似**：见 §2.3。
- **基准率**：训练集 1:1，与 genome-wide 约 1.69% 不同；报告 PPV 时已显式换算。
- **空间相关**：按染色体留出可防止跨染色体泄漏，但同一染色体上的窗口仍空间相关，bootstrap CI 可能偏窄。
- **单菌株/单组装**：结论不外推到其他菌株或近缘种。

---

## 7. 许可与引用

- 代码：MIT（见 `LICENSE`）。
- 数据：不随仓库分发；请遵守各来源条款。
- 建议引用：
  - OriDB：Nieduszynski CA, Knox Y, Donaldson AD (2012). OriDB: a DNA replication origin database. *Nucleic Acids Research*.
  - SGD：Saccharomyces Genome Database (Stanford University).
  - RefSeq：O'Leary NA et al. (2016). Reference sequence (RefSeq) database at NCBI. *Nucleic Acids Research*.
- 本项目仅用于统计与方法教学演示，不构成实验结论。

---

## 8. 离线使用

若无法联网：

1. 用 `data/sample/sample_dataset.csv`（45 个窗口，含坐标、序列、特征）运行 `python scripts/smoke_test.py`，验证代码可执行；
2. 正式分析需要按 §2 的 URL 手工下载三个文件，放入 `data/raw/`，文件名必须与 `ars/config.py` 一致：
   - `GCF_000146045.2_R64_genomic.fna.gz`
   - `SGD_features.tab`
   - `oridb_sc_ori.tsv`（可由 OriDB 的 `data_output.php` 表格另存为 TSV）
3. 然后执行 `python scripts/run_all.py --skip-download`。

示例数据只用于冒烟测试，**不得**用于报告正式 AUC 或结论。
