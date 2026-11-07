# 书中统计知识点 ↔ ARS 预测项目环节对照

> 项目：预测酿酒酵母（*S. cerevisiae*）复制起始点（ARS）的二分类小项目
> 教材：*Statistical Thinking for the 21st Century*（Russell A. Poldrack）
> 本文件只做「书—项目」的对应关系梳理，不包含项目代码与实验结果。

---

## 1. 书籍信息

| 项目 | 内容 |
| --- | --- |
| 书名 | *Statistical Thinking for the 21st Century* |
| 作者 | Russell A. Poldrack |
| 版权 | Copyright 2019 Russell A. Poldrack |
| 文件 | `StatsThinking21-2024.pdf`，5.6 MB |
| PDF 生成时间 | 2024-03-31（元数据 creationDate/modDate；creator 为 LaTeX via pandoc，producer 为 pdfTeX） |
| 页数 | 321 页（正文 p.9 起，References p.319–321） |
| 章节规模 | 前言 + 第 1–18 章 + 参考文献；PDF 目录共 157 条（20 条一级、137 条二级） |
| 性质 | 开放获取教材，源自作者在 Stanford 开设的 Psych 10 / Stats 60 课程（Preface, p.9） |
| 代码与数据 | 示例用 R 与 NHANES 数据集（如 Ch.5、Ch.7、Ch.14、Ch.17） |

**目录核对方式**：用 PyMuPDF 读取 PDF 元数据与内嵌书签目录（`get_toc()`，157 条），并对 p.9、p.238 等页做了印刷页码抽查，确认**下文页码均为 PDF 页码，且与书内印刷页码一致**。

**编号说明**：PDF 书签的一级条目不带章号，但正文页眉按 `Chapter 1 … Chapter 18` 编号，本文件按正文页眉编号。

---

## 2. 完整章节清单（以 PDF 实际目录为准）

| 章 | 标题（原文 / 中文） | 页码范围 | 二级目录（页码） |
| --- | --- | --- | --- |
| 前言 | Preface / 前言 | 9–12 | Why does this book exist? (9)；The golden age of data (10)；The importance of doing statistics (11)；An open source book (11)；Acknowledgements (12) |
| 1 | Introduction / 引言 | 13–22 | What is statistical thinking? (13)；Dealing with statistics anxiety (14)；What can statistics do for us? (14)；The big ideas of statistics (18)；Causality and statistics (20)；Learning objectives (21)；Suggested readings (22) |
| 2 | Working with data / 数据基础 | 23–32 | What are data? (23)；Discrete versus continuous measurements (25)；What makes a good measurement? (26)；Learning Objectives (29)；Suggested readings (29)；Appendix (29) |
| 3 | Summarizing data / 数据汇总 | 33–48 | Why summarize data? (33)；Summarizing data using tables (34)；Idealized representations of distributions (43)；Learning objectives (46)；Suggested readings (47) |
| 4 | Data Visualization / 数据可视化 | 49–64 | Anatomy of a plot (51)；Principles of good visualization (53)；Accommodating human limitations (58)；Correcting for other factors (61)；Learning objectives (62)；Suggested readings and videos (63) |
| 5 | Fitting models to data / 用模型拟合数据 | 65–96 | What is a model? (65)；Statistical modeling: An example (66)；What makes a model “good”? (75)；Can a model be too good? (77)；Summarizing data using the mean (79)；The mode (82)；Variability: How well does the mean fit the data? (83)；Using simulations to understand statistics (84)；Z-scores (85)；Learning objectives (94)；Appendix (94) |
| 6 | Probability / 概率 | 97–120 | What is probability? (97)；How do we determine probabilities? (99)；Probability distributions (106)；Conditional probability (108)；Computing conditional probabilities from data (110)；Independence (111)；Reversing a conditional probability: Bayes’ rule (113)；Learning from data (115)；Odds and odds ratios (116)；What do probabilities mean? (117)；Learning objectives (118)；Suggested readings (118)；Appendix (119) |
| 7 | Sampling / 抽样 | 121–128 | How do we sample? (121)；Sampling error (122)；Standard error of the mean (124)；The Central Limit Theorem (125)；Learning objectives (127)；Suggested readings (127) |
| 8 | Resampling and simulation / 重采样与模拟 | 129–138 | Monte Carlo simulation (129)；Randomness in statistics (130)；Generating random numbers (131)；Using Monte Carlo simulation (133)；Using simulation for statistics: The bootstrap (135)；Learning objectives (137)；Suggested readings (137) |
| 9 | Hypothesis testing / 假设检验 | 139–162 | Null Hypothesis Statistical Testing (NHST) (139)；Null hypothesis statistical testing: An example (140)；The process of null hypothesis testing (141)（含 9.3.5.3 Computing p-values using randomization，p.150）；NHST in a modern context: Multiple testing (158)；Learning objectives (161)；Suggested readings (161) |
| 10 | Quantifying effects and designing studies / 效应量化与研究设计 | 163–178 | Confidence intervals (163)；Effect sizes (170)；Statistical power (175)；Learning objectives (178)；Suggested readings (178) |
| 11 | Bayesian statistics / 贝叶斯统计 | 179–198 | Generative models (179)；Bayes’ theorem and inverse inference (180)；Doing Bayesian estimation (182)；Estimating posterior distributions (184)；Choosing a prior (189)；Bayesian hypothesis testing (191)；Learning objectives (196)；Suggested readings (196)；Appendix (196) |
| 12 | Modeling categorical relationships / 类别关系建模 | 199–208 | Example: Candy colors (199)；Pearson’s chi-squared test (200)；Contingency tables and the two-way test (202)；Standardized residuals (204)；Odds ratios (205)；Bayes factor (205)；Categorical analysis beyond the 2 X 2 table (206)；Beware of Simpson’s paradox (207)；Learning objectives (208)；Additional readings (208) |
| 13 | Modeling continuous relationships / 连续关系建模 | 209–222 | An example: Hate crimes and income inequality (209)；Is income inequality related to hate crimes? (210)；Covariance and correlation (210)；Correlation and causation (215)；Learning objectives (218)；Suggested readings (218)；Appendix (218) |
| 14 | The General Linear Model / 一般线性模型 | 223–246 | Linear regression (225)；Fitting more complex models (231)；Interactions between variables (234)；Beyond linear predictors and outcomes (238)；Criticizing our model and checking assumptions (238)；What does “predict” really mean? (240)；Learning objectives (242)；Suggested readings (243)；Appendix (243)。另：正文小节 14.6.1 Cross-validation (241) 未进入 PDF 目录，但为本项目关键参考 |
| 15 | Comparing means / 均值比较 | 247–264 | Testing the value of a single mean (247)；Comparing two means (249)；The t-test as a linear model (251)；Bayes factor for mean differences (253)；Comparing paired observations (253)；Comparing more than two means (257)；Learning objectives (261)；Appendix (262) |
| 16 | Multivariate statistics / 多元统计 | 265–294 | Multivariate data: An example (266)；Visualizing multivariate data (267)；Clustering (269)；Dimensionality reduction (279)；Learning objectives (290)；Suggested readings (293) |
| 17 | Practical statistical modeling / 实用统计建模 | 295–304 | The process of statistical modeling (295)（7 步流程，p.295；含 17.1.7 多重比较与 Tukey，p.302；17.1.8 混杂因素，p.302）；Getting help (304) |
| 18 | Doing reproducible research / 可复现研究 | 305–318 | How we think science should work (305)；How science (sometimes) actually works (306)；The reproducibility crisis in science (307)；Questionable research practices (312)；Doing reproducible research (314)；Doing reproducible data analysis (316)；Conclusion: Doing better science (316)；Learning objectives (317)；Suggested Readings (317) |
| — | References / 参考文献 | 319–321 | — |

---

## 3. 项目环节 ↔ 书中章节对照表

> “书中支持度”一列区分三种情况：**直接**（书中有对应方法）、**概念**（书中只有相关概念，方法需自行迁移）、**缺口**（书中未覆盖，需外部资料）。这一列建议在报告的方法部分保留，避免把“自己扩展的做法”写成“书中的方法”。

| # | 项目环节 | 对应章节与页码 | 书中核心概念 | 本项目落地方式 | 书中支持度 |
| --- | --- | --- | --- | --- | --- |
| 1 | 数据获取与抽样 | Ch.7 Sampling, p.121–127（7.1–7.4）；Ch.2 Working with data, p.23–29；Ch.17 17.1.2–17.1.3, p.296–297 | 总体与样本；代表性抽样与偏差（bias）；抽样误差与标准误（SEM = σ̂/√n）；变量类型（定性/定量、离散/连续）；数据准备 | ARS 正样本来自注释库、负样本从基因组背景区域抽取；在报告中给出**正/负样本总数、正负比（如 1:1、1:2）、16 条染色体各自的样本数**；说明“负样本不是总体的随机样本”，因此存在选择偏差 | 概念 |
| 2 | 负样本构造（两种方案对比） | Ch.6 6.4–6.9, p.108–117；Ch.12 12.1–12.5、12.8, p.199–207 | 条件概率与独立性；基准率（base rate）；2×2 列联表与卡方检验；标准化残差；优势比（odds ratio）；Simpson 悖论 | 两种负样本（如“随机基因间区” vs “按 GC/长度匹配的对照”）各建一份数据集；用 2×2 表与卡方检验比较两套负样本的特征分布差异；用 OR 量化某个特征在两套负样本中的富集程度；用 Simpson 悖论提醒“合并/分染色体统计可能反转结论” | 概念 |
| 3 | 特征描述 | Ch.3 Summarizing data, p.33–47；Ch.5 5.5–5.9, p.79–93；Ch.4 Data Visualization, p.49–63 | 汇总表；分布形态（直方图/密度）；均值、中位数、众数；变异程度；z-score 标准化；可视化原则 | 输出“按类别（正/负）的特征描述表”：n、均值、中位数、SD、缺失率；对连续特征做 z-score 标准化后再入模；画正负类特征的密度图/箱线图（对应书中“先看分布再建模”的顺序） | 直接 |
| 4 | 假设检验 | Ch.9 9.3, p.141–158（随机化 p 值见 9.3.5.3, p.150–154） | NHST 六步：参数→H0/H1→检验统计量→零分布→p 值→显著性判定；t 统计量；用随机重排构造零分布；Type I/II 错误与 α/β | ①对单个序列特征做正负类差异检验（t 检验或随机化检验）；②对模型 AUC 做**标签置换检验**（打乱正负标签 ≥1000 次，得到零分布，p = (#随机 AUC ≥ 观测 AUC + 1)/(n+1)），这正是书中 9.3.5.3 的思路 | 直接 |
| 5 | 模型与回归 | Ch.14 14.1–14.4, p.225–238；Ch.5 5.1–5.4, p.65–79；Ch.13 13.3, p.210–215 | 模型是对数据生成过程的简化理论；参数估计与拟合优度；RMSE；过拟合（8 阶多项式例，p.77–79）；协方差与相关；广义线性模型（仅 14.4 一段提及） | 主模型为逻辑回归（二分类），随机森林作对照；报告系数/OR 及方向；用相关矩阵筛掉高度冗余的序列特征；用训练集/交叉验证两列指标展示过拟合 | 概念（逻辑回归未展开） |
| 6 | 交叉验证 | Ch.14 14.6.1, p.241–242；Ch.5 5.4, p.77–79；Ch.17 17.1.4–17.1.7, p.297–302 | 留出法思想；k 折交叉验证；用 CV 估计“新数据上的表现”；模型设定与模型批判 | **按染色体分组**的 5 折交叉验证（同一染色体的样本只出现在同一折），对比普通随机 5 折；报告每折 AUC、均值±SD，以及 train 与 CV 的差距 | 概念（分组 CV 属缺口） |
| 7 | 模型评估指标 | Ch.6 6.7, p.113–115；Ch.9 9.3.6.1、9.3.7, p.156–158；Ch.12 12.3, p.202–204；Ch.10 10.2, p.170–174 | 敏感度/特异度与基准率（PSA 例：敏感度 80%、特异度 70%、基准率 5.8% → 阳性后验概率仅 16%）；混淆矩阵视角；Type I/II 错误；优势比作为一种效应量 | 报告混淆矩阵（TP/FP/FN/TN 与正样本占比）、敏感度、特异度、精确率（PPV）、AUC、PR-AUC；**用 PSA 例说明为什么正样本稀少时精确率会很低**，从而解释为何要同时看 PR-AUC | 概念（AUC/PR-AUC 属缺口） |
| 8 | 不确定性 | Ch.8 8.5, p.135–137；Ch.10 10.1, p.163–170；Ch.7 7.3, p.124–125 | bootstrap：有放回重采样估计统计量的抽样分布；置信区间与边际误差；重叠 CI 不能当检验用（p.170）；SEM 与样本量的关系 | 对测试集预测做 **1000 次按染色体整组重采样的 bootstrap**，取 2.5%/97.5% 分位得到 AUC 与 PR-AUC 的 95% CI；同时报告点估计、CI 宽度、以及不同折之间的波动 | 概念（分组 bootstrap 属扩展） |
| 9 | 多重比较 | Ch.9 9.4, p.158–161；Ch.17 17.1.7, p.299–302；Ch.18 18.4, p.312–313 | 族错误率（familywise error）；Bonferroni 校正（α/m）（书中例：10⁶ 次检验、α=.05 → 期望 5×10⁴ 个假阳性；书中给出的校正后 α=0.000005，原文数值）；Tukey 方法；p-hacking | 设比较次数 m（特征数 × 负样本方案数，例如 10×2=20），校正阈值 0.05/m=0.0025；画“原始 p 值 + 校正阈值线”的图；报告校正前后显著特征数；在局限中说明未用 FDR | 直接（FDR 属缺口） |
| 10 | 结论与局限 | Ch.18 18.1–18.6, p.305–317；Ch.17 17.1.8, p.302–303；Ch.2 2.3, p.26–29；Ch.14 14.5, p.238–240；Ch.13 13.4, p.215–218；Ch.1 1.4, p.20–21 | 可复现性危机与 QRP（HARKing、p-hacking）；预注册；可复现分析；混杂因素；测量的信度与效度；模型假设检查（残差、Q-Q 图）；相关不等于因果 | 写“分析前计划”（预注册式：先定指标与折数，再跑）；固定随机种子并记录环境版本；**报告所有尝试过的负样本方案**（不只报最好的）；讨论标签质量（实验验证 ARS vs 计算预测）与混杂（染色体长度/GC 含量）；明确本项目只做预测不做因果 | 直接 |

---

## 4. 必须在报告里讲清楚的 8 个核心概念

> 下面每个概念都给出“书中出处 + 本项目要用的具体数字或图”。**数据相关数字需由 pipeline 实际输出填入**（用 `⟨…⟩` 标注占位），不要凭印象写。

### 4.1 总体 vs 样本、抽样偏差（Ch.7, p.121–127）

- **要点**：统计推断的前提是样本能代表总体；非随机抽样会带来系统性偏差，且增大样本量也消不掉偏差。
- **本项目呈现**：
  - 表：`n_pos = ⟨…⟩`、`n_neg = ⟨…⟩`、正负比 `⟨1:k⟩`；16 条染色体的样本量分布（最小/最大/中位数）。
  - 图：每染色体正负样本数堆叠柱状图。
  - 一句话：负样本是从“可注释背景区”抽的，不是全基因组的随机样本 → 存在选择偏差，AUC 的含义因此是“在该抽样框下”的排序能力。

### 4.2 基准率、敏感度、特异度与 PPV（Ch.6 §6.7, p.113–115）

- **要点**：书中 PSA 例子给出可复用的数字链：敏感度 80%、特异度 70%、基准率 5.8% → 阳性后验概率只有 16%，优势比 2.62。说明“检验看起来不错”不等于“阳性就可信”。
- **本项目呈现**：
  - 图：测试集混淆矩阵（TP/FP/FN/TN）。
  - 数字：正样本占比 `⟨prevalence⟩`、敏感度 `⟨…⟩`、特异度 `⟨…⟩`、精确率 `⟨…⟩`。
  - **关键论证**：两套负样本改变了正负比，敏感度/特异度基本不变，但精确率会明显变化 → 这是“两种负样本构造对比”最有说服力的一段。

### 4.3 过拟合与交叉验证（Ch.5 §5.4, p.77–79；Ch.14 §14.6.1, p.241–242）

- **要点**：书中 8 阶多项式在训练数据 RMSE=0.7，换新数据变成 1.1；12 折 CV 的 R²=0.60，既高于新数据的 0.34，又远低于原始数据的 0.95。结论：训练集指标必然虚高，CV 才是“新数据表现”的估计。
- **本项目呈现**：
  - 表：`train AUC`、`CV AUC（均值±SD）`、`每折 AUC`、`train−CV 差距`。
  - 判据：差距 > 0.05 视为过拟合信号（阈值需在报告中写明是自己设定的）。
  - 图：每折 AUC 的折线/点图（5 折 × 2 种 CV 方案）。

### 4.4 Bootstrap 与置信区间（Ch.8 §8.5, p.135–137；Ch.10 §10.1, p.163–170）

- **要点**：bootstrap 用“有放回重采样”近似统计量的抽样分布，不依赖正态假设；CI 给出估计的不确定性；**重叠的 CI 不能直接当作显著性检验**（p.170 明确提醒）。
- **本项目呈现**：
  - 数字：`1000` 次重采样、2.5%/97.5% 分位；AUC = `⟨点估计⟩` [`⟨下限⟩`, `⟨上限⟩`]；PR-AUC 同样给出 CI。
  - 图：bootstrap 分布直方图 + 点估计与 CI 竖线。
  - 方法说明：重采样单位是**染色体**（整组有放回），因为同染色体样本非独立；这一点书中没有，需注明是自行扩展。

### 4.5 效应量，而不只是 p 值（Ch.10 §10.2, p.170–174）

- **要点**：显著性 ≠ 重要性；Cohen's d 的经验档位 0.2/0.5/0.8（Table 10.1, p.171）；优势比也可作为效应量（p.174）。
- **本项目呈现**：
  - 表：每个序列特征一行，给 `Cohen's d` 或 `OR`、95% CI、样本量。
  - 图：森林图（效应量 + CI），按效应量排序。
  - 报告句式：“GC 含量每增加 1 个 SD，ARS 的优势比为 `⟨OR⟩`（95% CI `⟨…⟩`），属于 `⟨可忽略/小/中/大⟩` 效应”。

### 4.6 多重比较与 Bonferroni（Ch.9 §9.4, p.158–161）

- **要点**：检验次数越多，族错误率越高；书中用 10⁶ 次检验演示 α=.05 下期望 5×10⁴ 个假阳性，Bonferroni 把 α 降到 0.000005。
- **本项目呈现**：
  - 数字：比较次数 `m = ⟨…⟩`；校正阈值 `0.05/m = ⟨…⟩`。
  - 图：原始 p 值条形图 + 校正阈值水平线；标注“校正前 k 个显著 / 校正后 j 个显著”。
  - 说明：本项目用 Bonferroni（保守），FDR 留作后续。

### 4.7 p 值的正确解释（Ch.9 §9.3.7, p.157–158）

- **要点**：书中直接列出三种常见误读：p 不是“重复实验得到同样结果的概率”，也不是“效应重要”的证据。
- **本项目呈现**：
  - 在结果段固定一句模板：“在 H0（正负样本无可区分信号）下，出现不低于观测 AUC 的概率为 `⟨p⟩`”，而不是“模型正确的概率是 1−p”。
  - 每个 p 值旁边同时给效应量/CI（见 4.5），避免只靠显著性下结论。

### 4.8 可复现性与可疑研究实践（Ch.18, p.305–317）

- **要点**：书中列出的 QRP 包括 HARKing、按 p 值挑变量、逐步加样本直到显著、随意剔除样本、变换数据直到显著；对策是预注册与可复现分析。
- **本项目呈现**：
  - 数字/记录：随机种子 `⟨seed⟩`；Python 与主要包版本表；`n` 次实验（含失败方案）全记录。
  - 图/表：“分析计划 vs 实际执行”对照表（计划改动的每一处都要写明原因）。
  - 明确写出：两种负样本方案**都报**，不因某个方案 AUC 低而删掉。

---

## 5. 本书不覆盖、可作为“局限与后续”的内容

以下内容已用全文关键词检索确认：书中**没有**专门章节或系统论述。写进报告时不要标注为“书中方法”。

| 缺口 | 书中实际情况 | 建议写入报告的表述（局限/后续） |
| --- | --- | --- |
| 逻辑回归 | 仅 Ch.14 §14.4（p.238）一句“可用广义线性模型处理二分类结果，后文再展开”；但全文检索无 logistic regression 章节，后文并未展开 | 逻辑回归的损失函数、正则化、系数/OR 解释需引用外部教材；本书只提供 GLM 与“模型—拟合—批判”的框架 |
| 随机森林 / 树模型 / 集成学习 | 未覆盖；Ch.16（p.266）只区分了 supervised 与 unsupervised learning | 随机森林作为“对照模型”属自行引入；报告中需外部引用，并说明它只用于检验“非线性/交互是否带来增益” |
| ROC / AUC / PR-AUC / 阈值选择 / 概率校准 | 未覆盖；全书只有敏感度、特异度、PPV、基准率与效应量 | AUC 的“随机正负样本对排序正确率”解释、PR-AUC 在不平衡数据下的优势、阈值选择与校准，均需外部引用；本书可用于解释**为什么单一阈值指标会误导** |
| 类别不平衡与负样本采样策略 | 未覆盖；只有基准率、优势比、列联表等概念 | “两种负样本构造”的动机（分布偏移、伪负样本）需自行论证；可用 Ch.6/Ch.12 的基准率与列联表逻辑支撑评估部分，但采样策略本身属外部知识 |
| 分组/分块交叉验证（按染色体） | 未覆盖；Ch.14 §14.6.1（p.241–242）的 CV 默认样本独立同分布 | 需说明同染色体 ARS 在序列组成与复制时间上相关，普通随机折会造成**信息泄漏**；分组 CV 是自行扩展，可用 Ch.7 的“代表性抽样/独立性”作为概念铺垫 |
| 序列特征工程（k-mer、GC 含量、位置/间距等） | 未覆盖；Ch.3/Ch.4 只给描述统计与可视化通则 | 特征构造与生物学解释需外部引用；本书只负责“如何描述与可视化这些特征” |
| FDR / Benjamini–Hochberg | 未覆盖；Ch.9 只有 Bonferroni，Ch.17 只有 Tukey | 明确本项目采用保守的 Bonferroni，FDR 控制列为后续工作 |
| 功效分析（power analysis） | 有（Ch.10 §10.3, p.175–178），但本项目未做事前功效计算 | 可在“后续”中补：给定预期效应量与 α，估算所需正/负样本量；书中 n=64/组、d=0.5、power=0.8 的例子可直接引用作模板 |
| 贝叶斯方法 / Bayes factor | 有（Ch.11, p.179–198；Ch.12 §12.6, p.205）但本项目未使用 | 可列为后续：用 Bayes factor 替代 p 值做模型比较；或对 AUC 给后验分布 |
| 聚类与降维 | 有（Ch.16, p.269–290）但本项目未使用 | 可列为后续：用特征聚类去冗余、用 PCA 检查特征共线性 |
| 因果推断 | Ch.1 §1.4（p.20–21）与 Ch.13 §13.4（p.215–218）明确讨论“相关≠因果” | 明确本项目是**预测**任务：AUC 高不代表某个序列特征“导致”ARS 活性；这也是报告结论部分必须写的一句 |
| 测量效度/信度 | Ch.2 §2.3（p.26–29）有框架（构念、信度、效度） | 可用于讨论标签质量：实验验证的 ARS 与计算预测的 ARS 属于不同“测量”，存在标签噪声；本项目未做标签一致性评估，列为局限 |
| 工具与数据差异 | 全书用 R + NHANES 示例数据 | 书中代码不能直接复用；报告需说明 Python 实现与原书示例的差异，避免读者误以为可直接对照复现 |

---

## 6. 引用与写作提醒

1. **页码格式**：统一写“Ch.X §X.Y, p.Z”，其中 p.Z 为 PDF 页码（已核对与印刷页码一致）。
2. **三种支持度要区分**：报告中凡属“书中直接方法”的（bootstrap、随机化 p 值、Bonferroni、CV、CI、效应量），可标注书中页码；凡属“自行扩展”的（分组 CV、分组 bootstrap、AUC/PR-AUC、逻辑回归细节、随机森林），必须写明“本书未覆盖，参考外部资料”，否则容易被评审认为误引。
3. **不要用书的结论替代本项目结论**：书的例子是 NHANES/R 场景，数字（如 RMSE 0.7→1.1、R² 0.95/0.34/0.60、PSA 16%、α=0.000005）只能作为**方法论类比**，不能当作本项目的实验结果。
4. **报告结构建议**（与第 3 节顺序一致）：数据与抽样 → 负样本构造 → 特征描述 → 单特征检验 → 模型 → 交叉验证 → 评估指标 → 不确定性 → 多重比较 → 可复现性与局限。这一顺序也对应 Ch.17 §17.1 的 7 步建模流程（p.295）。
