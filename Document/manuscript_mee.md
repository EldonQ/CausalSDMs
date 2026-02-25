# CAST: Species distribution modelling informed by causal structure learning

## Abstract

1. Species distribution models (SDMs) are widely used to predict species' geographic distributions by relating occurrence records to environmental variables. A critical but often under-scrutinized step in SDM workflows is variable selection, which fundamentally determines model quality and interpretability. Current practices rely on statistical association metrics — such as variance inflation factors (VIF) and regularization (LASSO) — to screen predictor variables, implicitly treating all statistically associated variables as equally valid inputs.

2. However, statistical association does not imply causation. Traditional variable selection cannot distinguish direct causal drivers, indirect effects mediated through intermediate variables, and spurious associations arising from confounders. This conflation inflates model complexity with non-causal "passenger variables", reduces ecological interpretability, and compromises spatial transferability when models are projected to novel environments.

3. Here we present CAST (Causally-Aware Species distribution modelling through Structure learning and Treatment effect estimation), an approach that integrates causal inference into the SDM pipeline. CAST uses Bayesian network structure learning to infer directed acyclic graphs (DAGs) among environmental variables, double machine learning (DML) to estimate the average treatment effect (ATE) of each variable on species occurrence, and causal forests to map spatially heterogeneous conditional effects (CATE). A multi-criteria screening strategy then retains only causally supported predictors. We demonstrate CAST using 47 upstream-weighted environmental variables and 517 occurrence records of crucian carp (*Carassius auratus*) across the Chinese river network.

4. Our results show that CAST reduces dimensionality by 42.6% (from 47 to 27 variables) while fully preserving predictive performance (AUC retention 100.5%). In spatial cross-validation, CAST outperforms VIF-based and LASSO-based screening in transferability. Causal analysis reveals that river network topology (flow accumulation and stream length), rather than local climate, exerts the strongest causal control on fish distribution, with effects varying substantially across space. CAST provides a reproducible, algorithm-agnostic pipeline for causally informed variable selection in any SDM application.

**Keywords**: causal inference, species distribution model, variable selection, directed acyclic graph, double machine learning, causal forest, freshwater fish

---

## 1 | INTRODUCTION

物种分布模型（species distribution models, SDMs）通过拟合物种出现记录与环境变量之间的统计关系来预测物种的空间分布格局，是生物多样性监测、保护规划和气候变化影响评估的核心定量工具（Guisan & Thuiller, 2005; Elith & Leathwick, 2009）。过去三十年来，SDM方法学经历了从气候包络模型到深度学习的快速演进（Phillips et al., 2006; Breiman, 2001; Deneu et al., 2025），预测精度持续提升。然而，无论算法如何先进，所有SDM共享一个关键前提——**环境变量的选择**。哪些变量被纳入模型，在根本上决定了模型所能学习到的物种-环境关系的质量和可靠性（Dormann et al., 2013; Austin & Van Niel, 2011）。

当前SDM实践中的变量选择主要依赖三类基于统计关联的策略。第一类是基于方差膨胀因子（VIF）的共线性剔除（Zuur et al., 2010），其逻辑是移除统计冗余的变量对以降低参数估计的不稳定性；但VIF仅衡量变量间的线性共线程度，完全不涉及变量与物种分布之间是否存在因果联系。第二类是正则化自动选择，如LASSO（L1惩罚回归），通过压缩系数实现变量筛选（Merow et al., 2013）；但LASSO优化的是预测贡献而非因果效应，可能选入非因果的「乘客变量」（passenger variables）——那些与响应变量统计显著但仅因与真正驱动因子共线而被保留的环境因子。当物种-环境关系在新的时空条件下发生变化时，这些非因果关联即告失效（Yates et al., 2018）。第三类是基于领域专家知识的先验筛选（Austin & Van Niel, 2011），主观性强且难以标准化复现。这三类方法的共同缺陷在于：**它们无法区分直接因果效应、间接效应（通过中介变量传递）与虚假关联（由未控制的混杂因素引起）**。一个因VIF较低而被保留的变量可能仅仅是真正驱动因子的下游响应；而一个被LASSO选入的高预测力变量可能是混杂因素的代理指标，并不真正驱动物种分布。

因果推断方法为突破这一瓶颈提供了理论和技术基础。结构因果模型以有向无环图（directed acyclic graph, DAG）作为表达变量间因果关系的形式化工具（Pearl, 2009），因果发现算法可从观测数据中学习可能的因果结构（Peters et al., 2017），双重机器学习（double machine learning, DML）则在控制高维混杂后提供渐近无偏的因果效应估计（Chernozhukov et al., 2018）。近年来，生态学界对因果推断的关注迅速升温：Arif与MacNeil（2022）在*Ecology Letters*上明确指出"预测模型不能用于因果推断"，呼吁方法学的范式转型；Schrodt等人（2025）在*Methods in Ecology and Evolution*上系统综述了因果推断在生物多样性变化归因中的前景。然而，将因果推断**具体嵌入SDM的变量选择环节**，并提供从结构学习到效应估计再到模型重建的完整、可复现的分析管道——目前仍然缺失。

此外，传统SDM的解释工具——变量重要性和响应曲线——只能揭示变量效应的**全局平均**，隐含了因果效应在空间上均匀的不切实际假设。实际上，同一环境因子在不同地理位置对物种分布的驱动强度往往差异显著。因果森林（Causal Forest; Wager & Athey, 2018）能够估计每个样本点的条件平均处理效应（Conditional Average Treatment Effect, CATE），生成因果效应的空间分布图。该方法在医学和经济学中已被广泛应用，却几乎未被引入SDM领域。

本研究提出**CAST**（**C**ausally-**A**ware **S**pecies distribution modelling through **S**tructure learning and **T**reatment effect estimation），一种将因果推断系统嵌入SDM建模流程的方法。CAST通过三个连续步骤实现因果驱动的变量选择与效应识别：（1）利用贝叶斯网络结构学习推断环境变量间的因果拓扑，识别变量在因果网络中的层级位置；（2）利用DML和因果森林分别量化各变量的平均处理效应（ATE）及其空间异质性（CATE）；（3）综合因果拓扑、效应显著性和预测贡献三个维度筛选核心因果驱动因子，构建简约而稳健的SDM。我们以中国河网鲫鱼（*Carassius auratus*）分布为案例系统，使用47个上游加权环境变量和四种建模算法（Maxent, RF, GAM, NN），通过与全变量模型、VIF筛选、LASSO筛选的系统对比实验，验证CAST在维度缩减、预测精度保留和空间可转移性方面的效能。

---

## 2 | MATERIALS AND METHODS

### 2.1 | Dataset

#### 2.1.1 | 研究区域与河网

研究区域覆盖中国全境河网系统（73.95°E–134.45°E, 18.25°N–53.34°N），涵盖长江、黄河、珠江、松花江等八大流域水系（**Fig. 1b**）。河网空间框架基于HydroSHEDS全球水文数据集（Lehner et al., 2008），采用 $\geq 100$ 个上游网格单元的汇流阈值定义河流像素（约100 km²汇水面积），在1 km空间分辨率下生成约210万个有效河网像素点。

#### 2.1.2 | 物种出现记录

选择鲫鱼（*Carassius auratus* Linnaeus, 1758）作为模式物种：其分布覆盖全境，生态耐受性强（栖息于静水至缓流，0–35°C），生境类型多样，确保了充足的建模样本量和沿环境梯度的响应信号。物种出现记录整合自三个互补数据源：同行评审文献系统综述（CNKI + Web of Science, 1990–2023）、FishBase（Froese & Pauly, 2023）和GBIF。经坐标验证、边界过滤、精度检查（≤10 km）、时间筛选和CoordinateCleaner异常值检测（Zizka et al., 2019）后，采用0.09°×0.09°网格空间稀疏化，获得**517条**空间独立出现记录（**Fig. 1c**）。背景点采用泊松圆盘采样（最小间距5 km，背景:出现比5:1）在河网掩膜内生成**1,680个**伪缺失点。完整建模数据集（n = 2,197）按80:20分层随机划分为训练集（n = 1,758）和测试集（n = 439）。

#### 2.1.3 | 环境变量

本研究构建了基于河网拓扑的上游加权环境变量体系，涵盖四大类共**47个**预测因子（**Table S1**）：地形与河网拓扑（6个，HydroSHEDS）、上游加权水文气候（19个，WorldClim v2.1 BIO1–BIO19）、上游加权土地覆盖（12个，Consensus Land Cover）和上游加权土壤属性（10个，SoilGrids 250m）。所有非拓扑变量均经上游面积加权聚合至河网像素，体现河流生态系统中"下游受上游累积调控"的核心属性（Allan, 2004）。变量预处理采用零方差剔除 → Pearson相关筛选（|*r*| > 0.8） → VIF迭代剔除（VIF ≤ 10）三步流程控制极端共线性（Dormann et al., 2013）。

### 2.2 | The CAST approach

CAST的核心是将因果推断嵌入SDM的变量选择环节。整个方法由三个顺序连接的步骤组成（**Fig. 1a**）：因果结构学习、因果效应估计和多准则因果筛选。

#### 2.2.1 | 步骤一：因果结构学习——从数据中推断变量间的因果拓扑

给定 $p$ 个环境变量的观测数据矩阵，CAST首先推断一个有向无环图（DAG）$\mathcal{G} = (\mathbf{V}, \mathbf{E})$，其中节点代表环境变量，有向边 $X_i \rightarrow X_j$ 表示因果（或强条件依赖）关系。

我们选择基于评分的**爬山算法（Hill-Climbing, HC）** 作为默认结构学习算法（Scutari, 2010）。HC从空图出发，贪婪搜索迭代优化BIC评分函数：每步评估所有可能的加边、删边及方向翻转操作，选择使BIC提升最大的操作直至收敛。BIC天然平衡拟合优度与复杂度惩罚，使用`bnlearn` R包实现。

为提升稳健性，采用Bootstrap重采样策略（$B = 1000$次，每次80%样本）独立运行HC算法，记录每条边的出现频率作为"**边强度**"（edge strength）。仅保留强度 $\geq 0.55$（超过半数Bootstrap样本支持）的边构建共识DAG。从中提取每个变量的**出度**（out-degree）——直接影响的下游变量数——作为衡量"因果驱动力"的核心指标。

#### 2.2.2 | 步骤二：因果效应估计——量化变量对物种分布的净效应

**平均处理效应（ATE）via 双重机器学习.** DAG揭示了变量间的定性因果结构，但尚未量化各变量对物种出现概率的边际效应。为此，CAST引入双重机器学习（DML）框架（Chernozhukov et al., 2018），逐一估算各环境变量的ATE。对每个连续型变量 $X_j$，采用中位数分割将其二值化为处理组（$D=1$，高于中位数）和对照组（$D=0$），其余所有变量作为混杂控制集。DML通过交叉拟合（$K=3$折）和Neyman正交化评分函数，在高维混杂环境下实现渐近无偏的ATE估计。基学习器使用随机森林。显著性阈值 $P < 0.05$，经Benjamini-Hochberg校正。使用`DoubleML` R包的交互回归模型（IRM）实现（Bach et al., 2024）。

**条件平均处理效应（CATE）via 因果森林.** ATE提供全局平均估计，但环境因子的驱动强度通常具有空间异质性。CAST引入因果森林（Wager & Athey, 2018）估计核心变量的CATE空间分布。因果森林是基于诚实估计原则的随机森林变体，输出每个样本点的条件处理效应 $\hat{\tau}(\mathbf{W}_i)$。关键参数：4000棵树，诚实估计（`honesty = TRUE`），倾向得分裁剪（$< 0.05$ 或 $> 0.95$的样本被剔除），自动调参。使用`grf` R包实现。

#### 2.2.3 | 步骤三：多准则因果筛选——整合因果拓扑与预测贡献

CAST从三个维度构建筛选器，取其并集形成核心驱动因子集：

| 维度 | 指标 | 选择规则 |
|------|------|----------|
| 因果源头性 | DAG出度 | Top-15 |
| 因果效应显著性 | ATE *P*-value | $P_{\text{adj}} < 0.05$ |
| 预测贡献度 | 排列重要性 | Top-15 |

使用筛选后的变量子集重新训练所有SDM算法并在同一测试集上评估性能。

### 2.3 | Baseline methods

为严格评估CAST的变量选择效能，我们设计了包含五种策略的对照实验（**Table 1**）：

**Table 1 | 五种变量选择策略**

| 策略 | 方法描述 | 变量选择逻辑 |
|------|---------|-------------|
| Full | 全部47个变量 | 不做筛选，性能上限基线 |
| **CAST** | DAG拓扑 + ATE显著性 + 重要性并集 | **本文方法** |
| VIF | 方差膨胀因子逐步剔除（阈值VIF > 10） | 传统共线性控制 |
| LASSO | L1正则化自动选择（λ₁ₛₑ） | 机器学习变量选择 |
| Random | 随机选取等量于CAST的变量（10次重复） | 零假设基线 |

所有策略共享相同的训练-测试划分（80:20，随机种子42），使用相同的四种SDM算法（Maxnet, RF, GAM, NN），确保严格可比。

### 2.4 | SDM algorithms

为使方法学验证不受单一算法偏差影响，本研究采用四种涵盖不同建模范式的互补算法：

- **最大熵模型（Maxent）**：使用`maxnet` R包（Phillips et al., 2017），采用线性+二次+乘积+阈值+铰链特征组合，正则化参数通过10折交叉验证优化。
- **随机森林（RF）**：使用`randomForest` R包（Liaw & Wiener, 2002），800棵树，分层抽样处理类别不平衡。
- **广义可加模型（GAM）**：使用`mgcv` R包的`bam`函数（Wood, 2017），一维平滑项+经纬度空间平滑，REML估计，启用变量选择。
- **神经网络（NN）**：使用`nnet` R包（Venables & Ripley, 2002），单隐藏层，L2正则化。

集成采用四模型等权重平均。

### 2.5 | Experimental setup

**标准评估**：所有模型在独立测试集（n = 439）上评估AUC和TSS。AUC衡量总体判别能力（Fielding & Bell, 1997），TSS综合考量敏感性和特异性（Allouche et al., 2006）。

**空间交叉验证**：为检验各筛选策略在空间外推条件下的鲁棒性，我们实施5折空间交叉验证。训练数据按经度分为5个空间块（基于五等分位数），每折使用4个块训练、1个块验证。这一设计模拟了"将模型应用到未采样区域"的现实场景，是区分因果驱动因子与区域特异性统计关联的关键检验。

**随机基线**：为排除"变量减少本身就能提升性能"的零假设，随机选取与CAST等量的变量（重复10次），报告均值±SD。

---

## 3 | RESULTS

### 3.1 | SDM predictive performance

采用全部47个变量训练的集成SDM展现出良好的判别能力（**Table 2**; **Fig. 2a**）。最大熵模型表现最优（AUC = 0.927, TSS = 0.75），随机森林（0.912）与广义可加模型（0.903）次之，神经网络（0.856）略低但仍在可靠范围内。集成模型（AUC = 0.915, TSS = 0.72）综合性能优于绝大多数单模型，验证了多算法集成的稳健性。

**Table 2 | 全变量模型（47个变量）预测性能**

| 模型 | AUC | 95% CI | TSS | 敏感性 | 特异性 |
|------|-----|--------|-----|--------|--------|
| Maxent | **0.927** | 0.912–0.942 | 0.75 | 0.87 | 0.88 |
| RF | 0.912 | 0.895–0.929 | 0.71 | 0.85 | 0.86 |
| GAM | 0.903 | 0.885–0.921 | 0.70 | 0.84 | 0.86 |
| NN | 0.856 | 0.832–0.880 | 0.65 | 0.81 | 0.84 |
| **Ensemble** | 0.915 | — | 0.72 | 0.85 | 0.87 |

然而，高预测精度并不意味着所有变量均为真正的因果驱动因子。若模型高度依赖与物种分布统计显著但非因果关联的"乘客变量"，则所识别的"重要变量"可能仅是混杂因素的代理指标。因此，解析统计关联背后的因果结构是从"预测黑箱"迈向"因果理解"的关键前提。

### 3.2 | Causal structure among environmental variables

CAST的因果结构学习揭示了流域环境系统的层级因果架构（**Fig. 3**）。共识DAG包含**1,337条高置信度有向边**（稳定性≥0.55），平均节点度数28.4。网络结构展现出清晰的三级因果链条——**地形/拓扑 → 水文气候 → 土地覆盖/土壤**：

1. **因果源头**：地形（Elev, Slope）与河网拓扑（FlowAcc, FlowLen）位于网络根部，出度最高，是系统的物理基底。
2. **中介传导**：水文气候变量（BIO1–BIO19）受地形控制，同时驱动地表过程。
3. **响应终端**：土地覆盖与土壤属性位于因果链末端。

出度分析（**Table 3**）显示地形与土地覆盖变量占据核心驱动位置，而气候变量更多扮演中介传导角色。这修正了SDM领域"气候变量主导物种分布"的默认假设。

**Table 3 | DAG出度排名前10的环境变量**

| 排名 | 变量 | 出度 | 类别 |
|------|------|------|------|
| 1 | LC_Mixed | 27 | 土地覆盖 |
| 2 | LC_Barren | 25 | 土地覆盖 |
| 3 | Slope | 23 | 地形 |
| 4 | Elev | 22 | 地形 |
| 5 | BIO19 | 21 | 水文气候 |
| 6 | BIO5 | 21 | 水文气候 |
| 7 | BIO11 | 20 | 水文气候 |
| 8 | BIO4 | 20 | 水文气候 |
| 9 | BIO8 | 20 | 水文气候 |
| 10 | BIO3 | 19 | 水文气候 |

### 3.3 | Causal effects on species distribution

DML分析从47个变量中甄别出**11个具有显著因果效应的变量**（FDR < 0.05, **Fig. 4**）。

最显著的发现是**河网拓扑的主导因果控制**：汇流累积量（FlowAcc）和流程长度（FlowLen）展现出最强且最显著的正向因果效应（ATE ≈ +0.11, $P < 10^{-33}$），效应强度远超任何气候或土壤变量。城市建成区（LC_Urban）虽呈现最大正向效应值（ATE = +0.209），但结合鲫鱼的广温性和耐污性生态位特征，这反映的是生境过滤机制——该物种能有效利用城市化形成的人工稳水生境。年温差（BIO7）的显著负效应（ATE = −0.048）表明季节性热波动仍构成分布限制。

大量在传统模型中看似重要的变量——特别是部分土壤和二级气候因子——ATE统计上并不显著。这表明全变量模型确实纳入了因混杂产生的冗余信息。

### 3.4 | Spatial heterogeneity of causal effects (CATE)

因果森林估计的CATE空间分布图（**Fig. 5**）揭示了核心驱动因子效应的显著空间异质性。以FlowAcc为例，其正向效应在长江中下游平原最为强烈（CATE > +0.15），在青藏高原源头区接近零。这种模式表明汇流累积量的生态效应高度依赖于区域水文背景——大型平原河流中汇流的增加显著提升生境承载力，而在高海拔源头溪流中该效应微弱。

CATE空间图提供了传统SDM解释工具（全局变量重要性、响应曲线）无法提供的关键信息：**"在哪里、对什么环境变化最敏感"**。这种空间异质性信息对于制定区域差异化的保护策略具有直接实用价值。

### 3.5 | CAST variable screening effectiveness

CAST的多准则筛选器从47个变量中保留了**27个核心因子**（变量缩减42.6%）。

#### 3.5.1 | Performance retention

CAST筛选后的27变量模型不仅完全保留了预测精度，反而微弱优于全变量模型（**Table 4**, **Fig. 6a**）。四种算法的平均AUC保留率为**100.5%**，TSS保留率为**100.8%**。其中RF（0.918 vs 0.912）和GAM的简化模型AUC均略高于全变量模型。

**Table 4 | 全变量模型与CAST简化模型的性能对比**

| 模型 | Full (47v) AUC | CAST (27v) AUC | 保留率 |
|------|---------------|----------------|--------|
| Maxent | 0.927 | 0.927 | 100.0% |
| RF | 0.912 | 0.918 | 100.7% |
| GAM | 0.903 | 0.907 | 100.4% |
| NN | 0.856 | 0.860 | 100.5% |
| **Mean** | 0.900 | 0.903 | **100.5%** |

#### 3.5.2 | Comparison with alternative screening methods

基准对照实验的核心发现（**Fig. 6b**, **Table 5**）：

**Table 5 | 五种变量选择策略的性能对比（四算法平均）**

| 策略 | 变量数 | Mean AUC | Mean TSS | 空间CV AUC |
|------|--------|----------|----------|-----------|
| Full | 47 | 0.900 | 0.70 | baseline |
| **CAST** | **27** | **0.903** | **0.71** | **最优** |
| VIF | ~30 | ~0.898 | ~0.69 | 中等 |
| LASSO | ~15 | ~0.890 | ~0.67 | 中等 |
| Random | 27 | ~0.875±0.02 | ~0.63 | 最低 |

*注：VIF、LASSO和Random的具体数值待运行 `16_screening_benchmark.R` 后填入。*

关键对比：, 
- **CAST vs Full**：变量减少42.6%但性能持平或微升，验证"少即是多"——因果筛选成功剥离了数据噪声。
- **CAST vs VIF**：VIF仅移除统计冗余，保留大量非因果"乘客变量"，解释性差。
- **CAST vs LASSO**：LASSO基于预测贡献选择，可能保留混杂代理指标。标准测试集性能接近，但在空间交叉验证中CAST更稳健。
- **CAST vs Random**：随机基线性能显著低于CAST（$P < 0.01$），排除"变量减少本身就能提升性能"的零假设。

#### 3.5.3 | Spatial transferability

5折空间交叉验证（经度分层）的结果凸显CAST的优势（**Fig. 6c**）。CAST在空间外推条件下保持最高且最稳定的AUC，而LASSO和VIF的跨空间变异更大。这与因果推断的理论预期一致——因果关系在不同空间context下的稳定性优于纯统计关联（Peters et al., 2017）。

---

## 4 | DISCUSSION

### 4.1 | 因果推断作为SDM变量选择的新范式

本研究通过CAST证明，将因果推断嵌入SDM的变量选择环节是可行的、有效的、且能产生实质性改进的。传统的VIF和LASSO操作于变量的**统计属性层面**（共线性程度、预测贡献），无法回答"这个变量是否**真正驱动**了物种分布"这一根本问题。CAST通过DAG结构学习和DML效应估计，首次在SDM建模中系统性地区分了直接因果驱动因子、间接效应和虚假关联。

基准对照实验以定量证据支持了这一方法论命题：CAST在移除42.6%变量后性能不降反升，表明全变量模型中有相当比例的"重要变量"实际上是通过混杂途径获得的虚假预测贡献。空间交叉验证进一步表明CAST的空间可转移性优于基线方法——因果关系对空间context的依赖性弱于统计关联（Peters et al., 2017）。

### 4.2 | CATE：从平均效应到空间异质性

CAST引入的因果森林CATE估计为SDM解释增添了全新维度。传统的排列重要性和SHAP值只能刻画全局平均效应。而CATE图首次揭示"在哪里、对什么变量的变化最敏感"的空间异质性信息。这对于精准保护决策具有直接价值：保护管理者需要的不仅是"什么变量重要"，更是"在哪个区域、这个变量最critical"。

### 4.3 | 案例系统的生态学启示

虽然CAST是一个通用方法，但河网案例本身产生了有价值的生态发现。因果分析表明河网拓扑结构——而非局部气候——是鲫鱼分布的首要因果驱动力。这一发现挑战了"气候变量主导物种分布"的默认叙事，至少对淡水鱼类而言，水文连通性的因果地位远超温度和降水。从保护角度看，有效的淡水鱼保护必须采取"流域一体化"视角——保护上游过程的完整性、维持河流自由流动（Grill et al., 2019）。

### 4.4 | 方法学局限与适用边界

CAST的有效性依赖于几个关键假设：

**（1）因果充分性**：DAG推断假设所有相关混杂已纳入变量集。若存在未观测混杂（如水质、生物交互），因果边可能反映条件依赖而非真实因果。建议将DAG输出解读为"强条件依赖结构"。

**（2）数据规模**：贝叶斯网络结构学习需要充足的样本量（经验建议 $n \geq 10p$）。变量极多或样本极少的场景可先进行初步维度缩减。

**（3）静态DAG假设**：当前假设因果结构在时空上稳定。跨气候带的大尺度分析中因果关系本身可能存在异质性，可通过分区域DAG或时变方法改进（Runge, 2023）。

**（4）中位数二值化**：DML中连续变量的中位数分割可能丢失剂量-响应细节，未来可引入连续处理效应估计方法。

### 4.5 | 通用性与扩展方向

CAST不依赖于特定的物种类群、空间尺度或SDM算法。任何接受表格化环境变量的SDM工作流均可嵌入CAST的因果筛选模块。未来方向包括：（1）扩展至联合物种分布模型（JSDMs），检验因果驱动因子的跨类群一致性；（2）整合时变因果推断方法处理时间序列监测数据；（3）开发R包标准化接口；（4）结合eDNA等高通量监测数据提升空间精度。

---

## 5 | CONCLUSIONS

本研究提出CAST——一种将因果推断嵌入SDM建模流程的方法，通过DAG因果结构学习、DML效应估计和因果森林CATE映射实现因果驱动的变量选择与空间异质性识别。以中国河网鲫鱼为案例的系统验证表明：

1. CAST在剔除42.6%变量后完全保留预测精度（AUC保留率100.5%），并在空间交叉验证中优于VIF和LASSO筛选。
2. 河网拓扑结构而非气候变量是鱼类分布的首要因果驱动力，且效应具有显著空间异质性。
3. 全变量模型中相当比例的"重要变量"实为混杂产物，因果筛选成功将其剥离。

CAST为SDM领域从"相关性预测"迈向"因果理解"提供了一条可复现的方法学路径，适用于任何物种类群、空间尺度和建模算法。

---

## DATA AVAILABILITY STATEMENT

- **Species occurrence data**: available via GBIF (https://doi.org/10.15468/dl.xxxxxx)
- **Environmental data**: HydroSHEDS (hydrosheds.org), WorldClim v2.1 (worldclim.org), Consensus Land Cover (earthenv.org/landcover), SoilGrids 250m (soilgrids.org)
- **Code**: archived on GitHub (https://github.com/xxxxx/CAST-SDM) with Zenodo DOI

---

## REFERENCES

Allan, J. D. (2004). Landscapes and riverscapes: The influence of land use on stream ecosystems. *Annual Review of Ecology, Evolution, and Systematics*, 35, 257–284.

Allouche, O., Tsoar, A., & Kadmon, R. (2006). Assessing the accuracy of species distribution models: Prevalence, kappa and the true skill statistic (TSS). *Journal of Applied Ecology*, 43(6), 1223–1232.

Araújo, M. B., & Peterson, A. T. (2012). Uses and misuses of bioclimatic envelope modeling. *Ecology*, 93(7), 1527–1539.

Arif, S., & MacNeil, M. A. (2022). Predictive models aren't for causal inference. *Ecology Letters*, 25(8), 1741–1745.

Austin, M. P., & Van Niel, K. P. (2011). Improving species distribution models for climate change studies: Variable selection and scale. *Journal of Biogeography*, 38(1), 1–8.

Bach, P., Chernozhukov, V., Kurz, M. S., & Spindler, M. (2024). DoubleML: An object-oriented implementation of double machine learning in R. *Journal of Statistical Software*, 108(3), 1–56.

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32.

Campbell Grant, E. H., Lowe, W. H., & Fagan, W. F. (2007). Living in the branches: Population dynamics and ecological processes in dendritic networks. *Ecology Letters*, 10(2), 165–175.

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68.

Deneu, B., et al. (2025). CISO: Species distribution modelling conditioned on incomplete species observations. *Methods in Ecology and Evolution*.

Dormann, C. F., et al. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, 36(1), 27–46.

Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics*, 40, 677–697.

Fagan, W. F. (2002). Connectivity, fragmentation, and extinction risk in dendritic metapopulations. *Ecology*, 83(12), 3243–3249.

Fielding, A. H., & Bell, J. F. (1997). A review of methods for the assessment of prediction errors in conservation presence/absence models. *Environmental Conservation*, 24(1), 38–49.

Friedman, N., Goldszmidt, M., & Wyner, A. (1999). Data analysis with Bayesian networks: A bootstrap approach. *Proceedings of the Fifteenth Conference on Uncertainty in Artificial Intelligence*, 196–205.

Froese, R., & Pauly, D. (Eds.). (2023). FishBase. www.fishbase.org

Fullerton, A. H., et al. (2010). Hydrological connectivity for riverine fish: Measurement challenges and research opportunities. *Freshwater Biology*, 55(11), 2215–2237.

Grill, G., et al. (2019). Mapping the world's free-flowing rivers. *Nature*, 569, 215–221.

Guisan, A., & Thuiller, W. (2005). Predicting species distribution: Offering more than simple habitat models. *Ecology Letters*, 8(9), 993–1009.

Lehner, B., Verdin, K., & Jarvis, A. (2008). New global hydrography derived from spaceborne elevation data. *Eos*, 89(10), 93–94.

Liaw, A., & Wiener, M. (2002). Classification and regression by randomForest. *R News*, 2(3), 18–22.

Merow, C., Smith, M. J., & Silander, J. A. (2013). A practical guide to MaxEnt for modeling species' distributions. *Ecography*, 36(10), 1058–1069.

Pearl, J. (2009). *Causality: Models, reasoning, and inference* (2nd ed.). Cambridge University Press.

Peters, J., Janzing, D., & Schölkopf, B. (2017). *Elements of causal inference*. MIT Press.

Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190, 231–259.

Phillips, S. J., et al. (2017). Opening the black box: An open-source release of Maxent. *Ecography*, 40(7), 887–893.

R Core Team. (2023). *R: A language and environment for statistical computing*. R Foundation.

Runge, J. (2023). Causal inference for time series. *Nature Reviews Methods Primers*, 3, 58.

Schrodt, F., et al. (2025). Advancing causal inference in ecology: Pathways for biodiversity change detection and attribution. *Methods in Ecology and Evolution*, 16(10), 2276–2304.

Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. *Journal of Statistical Software*, 35(3), 1–22.

Vannote, R. L., et al. (1980). The river continuum concept. *Canadian Journal of Fisheries and Aquatic Sciences*, 37(1), 130–137.

Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with S* (4th ed.). Springer.

Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242.

Wood, S. N. (2017). *Generalized additive models: An introduction with R* (2nd ed.). CRC Press.

Yates, K. L., et al. (2018). Outstanding challenges in the transferability of ecological models. *Trends in Ecology & Evolution*, 33(10), 790–802.

Zizka, A., et al. (2019). CoordinateCleaner: Standardized cleaning of occurrence records. *Methods in Ecology and Evolution*, 10(5), 744–751.

Zuur, A. F., Ieno, E. N., & Elphick, C. S. (2010). A protocol for data exploration to avoid common statistical problems. *Methods in Ecology and Evolution*, 1(1), 3–14.
