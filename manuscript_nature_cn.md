# 因果推断揭示气候变化下机制驱动的淡水鱼类分布：基于河网视角

## 摘要

淡水生物多样性面临着气候变化日益严峻的威胁，然而现有的预测框架主要基于相关性，无法区分因果机制与虚假关联。本研究将因果发现与物种分布模型（SDM）相结合，在中国全境河网尺度上，利用 1 公里分辨率的河流专用环境图层进行了研究。以代表性淡水鱼类鲫鱼（*Carassius auratus*）为例，我们结合了四种 SDM 算法（Maxent、随机森林、GAM、神经网络）与基于约束和基于分数的因果结构学习算法，以绘制稳定的环境依赖关系。模型在独立测试数据上实现了稳健的判别能力（AUC 0.81–0.91）。因果网络分析识别了水文汇流、季节性水文气候变异性和地形梯度之间稳定的有向关系，并将其与单纯的相关性区分开来。通过因果森林（Causal Forests）估算的条件平均处理效应（CATE）揭示了空间异质的干预窗口，指出了栖息地适宜性对环境变化最敏感的区域。为了解决外推风险，我们仅使用“情景可用”变量重新训练模型，并在四种 SSP 路径下进行了未来预测（2041–2060）。预测的变化表明栖息地适宜性呈持续下降趋势，且结构不确定性（跨模型方差）超过了情景不确定性。我们的机制一致性框架使得能够优先考虑对气候敏感的河流路段进行适应性监测和保护，推动淡水生物多样性预测从基于相关性向因果知情转变。

**关键词：** 因果发现，物种分布模型，河网，淡水鱼类，气候变化，条件处理效应

## 1. 引言

淡水生态系统在不到地球表面 1% 的区域内蕴藏着不成比例的生物多样性，维持着超过 10% 的已描述物种和约三分之一的脊椎动物多样性（Dudgeon 等，2006；WWF，2020）。然而，受流量调节、污染、外来物种入侵、过度开发和气候快速变暖的综合影响，河流和湖泊中的种群数量下降和范围收缩速度快于陆地或海洋系统（Olden 等，2010；He 和 Silliman，2019；Reid 等，2019）。中国拥有世界上一些最多样化且受人类改造最严重的流域，在这里，气候变化与土地利用集约化及水利基础设施相互作用，重塑了淡水鱼类的栖息地条件。了解这些多重压力如何通过河网传播从而影响物种分布，不仅是一项生态挑战，也是紧迫的保护优先事项。

物种分布模型（SDMs）已成为绘制栖息地适宜性地图、预测气候驱动的范围转移以及为保护规划提供信息的主要量化框架（Guisan 和 Zimmermann，2000；Elith 和 Leathwick，2009；Guisan 等，2013；Peterson 等，2011）。相关性 SDM，包括 Maxent、随机森林和广义加性模型，可以捕捉物种出现与环境预测因子之间复杂的非线性关联，并且现在通常被组合成集合模型以提高稳健性（Araújo 和 New，2007；Thuiller 等，2009；Hao 等，2020；Valavi 等，2021）。然而，这些方法通常在纯粹的关联范式下运作。它们擅长识别共存的空间模式，但在揭示哪些环境变化直接主导分布转移方面提供的机制见解有限，并且在转移到新的环境条件（如未来气候）时往往表现不佳（Pearson 和 Dawson，2003；Yates 等，2018）。

这种局限性在河流生态系统中尤为突出，因为定向连通性、上游-下游累积效应和流域整合创造了标准网格化气候产品无法体现的强大空间依赖性（Poff 等，1997；Domisch 等，2015；Linke 等，2019）。对于受限于树状河网的水生生物而言，暴露度不仅由当地的气温或降水定义，还由上游土地利用、流量改变和地貌的累积遗产决定。然而，大多数淡水 SDM 应用仍然通过将河流线叠加在粗糙的陆地气候网格上来近似暴露度，掩盖了关键的因果路径，例如海拔梯度对温度的调节或聚集的农业径流导致的营养富集。同时，预测中的不确定性主要由算法和变量集之间的结构差异主导，但很少相对于情景不确定性进行明确分解（Lawler 等，2006；Thuiller，2004；Thuiller 等，2019）。

因果发现和基于机器学习的效应估算的最新进展为超越相关性、迈向具有机制意识的生物多样性预测提供了一条路径（Pearl，2000；Peters 等，2017；Ma 等，2020；Su 等，2021）。图因果模型和结构学习算法可以识别气候、水文、土地利用和地形之间稳定的依赖模式（Kalisch 和 Bühlmann，2007；Spirtes 等，2000），而双重机器学习和因果森林使得在高维设置下估算平均和异质处理效应成为可能（Chernozhukov 等，2018；Athey 和 Imbens，2016；Wager 和 Athey，2018）。在本研究中，我们将这些因果工具与多算法 SDM 集合相结合，开发了一个针对中国河网淡水鱼类分布的综合因果-预测框架，并以鲫鱼（*Carassius auratus*）为代表物种。具体而言，我们构建了 1 公里分辨率的河网一致性环境暴露图层，学习了连接水文、气候、地形和土地覆盖的共识因果网络，估算了关键驱动因素对出现率的平均和空间异质处理效应，并分解了本世纪中叶气候预测中的结构不确定性与情景不确定性。这些组件共同为优先考虑对气候敏感的河流路段进行适应性监测和保护提供了机制一致的基础。本研究的总体技术路线如图 1 所示。

**[插入图 1：研究技术路线图，展示数据组装、SDM 集合建模、因果发现及未来预测工作流程的整合。]**

## 2. 材料与方法

### 2.1 研究区域与物种出现数据

我们划定了研究区域，涵盖整个中国大陆河网系统（73.95°E–134.45°E, 18.25°N–53.34°N）。这一地理范围捕捉了极端的环境多样性，从东南部的亚热带季风区到西北部的干旱大陆性气候以及青藏高原的高山环境。该区域包括长江、黄河、珠江和松花江等主要流域，这些流域在温度、降水和人类足迹方面跨越了强烈的梯度。我们使用 ≥100 个单元格（约 100 km² 上游汇水面积）的汇流阈值定义河网像素，在 1 公里空间分辨率下产生了约 210 万个河流像素。以这种分辨率代表完整的河网，使我们能够同时捕捉到大型干流河道和可能作为淡水鱼类气候避难所或扩散走廊的小型支流。所有空间数据均投影到 Albers 圆锥等面积坐标系，以确保面积计算的准确性。

鲫鱼的出现记录汇编自三个互补来源以最大化空间覆盖：（i）从系统综述中提取的同行评审文献；（ii）FishBase；以及（iii）全球生物多样性信息网络（GBIF）。我们实施了严格的质量控制流程，包括坐标验证、国家边界空间过滤、精度检查和去重。为了减少机会性数据集中常见的采样偏差，我们应用了空间稀疏化，保留每个 0.09°（约 10 公里）网格单元中的一条记录。最终数据集包含 n=517 条空间独立的出现记录。这些出现记录在河网中的空间分布如图 2 所示。鲫鱼是一种体型较小的广温性鲤科鱼类，栖息于静水和流速缓慢的流水生境，包括水库、池塘和低地泛滥平原，并以其对富营养化和中度缺氧的耐受性而闻名。这些特征，加上其与人类改造水域的密切联系，使其成为研究气候和土地利用驱动的环境变化如何通过河网传播的易处理模式物种。为了近似可用的环境空间，我们使用限制在河网掩膜内的泊松圆盘采样生成了代表“可用”水生栖息地的背景点。我们强制执行 5 公里的最小点间距和 5:1 的背景与存在点比例，产生了 1680 个背景点，确保了可用环境的平衡代表性，同时限制了在密集采样区域的聚类（Wisz 等，2008；Comte 和 Olden，2021）。

**[插入图 2：展示 1 公里河网及鲫鱼出现记录分布的研究区域地图。]**

### 2.2 环境变量组装

我们从四个不同领域组装了环境预测因子，明确考虑了河网拓扑结构和上游-下游连通性。首先，我们计算了水文网络拓扑指标，包括汇流累积量、流路长度和河流等级，以量化流域尺度的连通性及在网络中的位置。其次，利用 EarthEnv-Streams 数据集，我们提取了按上游集水面积加权的月平均温度和降水，以及季度季节性指标。这种上游加权确保了环境值反映了流向某一点的整个流域的综合状况，而不仅仅是局部状况。

第三，我们从 SRTM 数字高程模型推导了地形梯度，将海拔、坡度和起伏等指标聚合到 1 公里河流像素分辨率。第四，我们使用 Consensus Land Cover 数据集和 SoilGrids250m 计算了上游加权的土地覆盖和土壤分数。这提供了代表上游流域内不同土地利用类型（如森林、城市、农业）比例和土壤属性（如粘土含量、pH 值）的变量。选择水文和气候变量是为了代表长期平均条件和年内变异性，这些条件预计会影响温带淡水鱼类的新陈代谢、繁殖和越冬生存，而土地覆盖和土壤属性则捕捉了可能改变栖息地质量和营养状况的人类改造和基质特性（Filipe 等，2013；Olden 等，2010）。

从最初的 100 多个候选变量中，我们通过零方差剔除、成对相关性筛选（|r| > 0.8）和阈值为 VIF ≤ 10 的迭代方差膨胀因子（VIF）剔除来减少多重共线性。这一严格的选择过程产生了 47 个用于最终建模的独立预测因子，如表 1 所总结。这一选择遵循 SDM 变量选择的最佳实践建议，以避免参数估计不稳定和预测因子重要性膨胀（Dormann 等，2013；Guisan 和 Zimmermann，2000）。

**[插入表 1：模型中使用的选定环境变量（水文、气候、地形、土地覆盖）汇总。]**

### 2.3 物种分布建模（SDM）

我们采用多算法集合方法来捕捉物种与其环境之间复杂的非线性关系。组合数据集（n=2016）被划分为分层的训练集（80%）和测试集（20%），保留了主要环境梯度上存在点和背景点的流行率。我们利用了四种涵盖从参数模型到机器学习模型连续体的不同算法：
1.  **最大熵（Maxent）**：使用 `maxnet` R 包建模，具有灵活的特征转换（线性、二次、乘积、阈值和铰链），以近似物种的潜在地理分布（Phillips 等，2006；Elith 等，2011；Merow 等，2013）。
2.  **随机森林（RF）**：由 500 棵回归树组成的集合，因其处理高维交互、非线性响应和预测因子间共线性的能力而被选中，同时对过拟合具有相对稳健性（Breiman，2001）。
3.  **广义加性模型（GAM）**：拟合带有惩罚薄板回归样条和空间张量积平滑项，以模拟非线性响应，同时保持平滑效应的部分可解释性（Hastie 和 Tibshirani，1986；Elith 和 Leathwick，2009）。
4.  **神经网络（NN）**：具有权重衰减的单隐藏层前馈网络，能够学习更受限模型难以捕捉的复杂模式和交互。

这些算法共同提供了关于物种-环境关系的互补视角：GAM 强调平滑、可解释的响应，RF 和 NN 捕捉复杂的交互作用，而 Maxent 非常适合具有灵活特征类的存在-背景数据（Elith 等，2011；Valavi 等，2021）。我们不是选择单一的“最佳”模型，而是解释跨算法的模式，并利用它们的联合行为来量化结构不确定性（Araújo 和 New，2007；Thuiller 等，2009）。

模型在独立测试集上使用受试者工作特征曲线下面积（AUC）、真实技巧统计量（TSS）、敏感性、特异性和连续 Boyce 指数进行评估，以确保稳健的性能评估。这些指标总结了跨阈值和流行率的判别能力，但高 AUC 或 TSS 值并不能自动保证可转移性或生态现实性，特别是在新的气候条件下（Araújo 和 Guisan，2006；Rodríguez-Rey 等，2019；Yates 等，2018）。因此，我们将它们与因果分析和不确定性划分结合使用，以避免过度解读任何单一算法或指标。

### 2.4 因果推断框架

为了超越相关性，我们将因果发现和效应估算整合到我们的工作流程中。
1.  **因果结构学习**：我们使用 PC 算法（基于约束）和爬山算法（基于分数）推断有向无环图（DAGs），以恢复气候、水文、地形、土地覆盖和物种出现之间稳定的依赖模式（Kalisch 和 Bühlmann，2007；Spirtes 等，2000；Peters 等，2017）。对于数据的每个 bootstrap 样本，我们估计一个 DAG，然后总结 300 次重复中的边缘稳定性，仅保留稳定性频率 ≥ 0.55 的边缘以构建共识因果网络。这种基于稳定性的过滤减少了采样噪声和高维设置中可能出现的微弱条件独立性的影响。
2.  **平均处理效应（ATE）**：我们使用双重机器学习（DML）来估算关键变量对物种存在的因果效应。DML 通过使用单独的机器学习模型分别预测结果和处理，然后进行正交化以去除一阶估计误差，从而允许在存在高维干扰参数（混杂因素）的情况下对因果参数进行近似无偏估计（Chernozhukov 等，2018；Ma 等，2020）。我们专注于可解释的处理变量，如城市土地覆盖比例、降水和上游土壤属性，这些对应于合理的管理手段。
3.  **条件平均处理效应（CATE）**：我们使用因果森林估算空间异质处理效应，该方法调整随机森林以估算协变量空间中因果效应的异质性（Athey 和 Imbens，2016；Wager 和 Athey，2018；Su 等，2021）。对于每个处理变量，CATE 表面揭示了该变量的边际变化将在何处产生最大（正向或负向）的出现概率变化，为空间针对性干预提供了依据。

### 2.5 未来气候预测与不确定性分解

我们在四种 CMIP6 共享社会经济路径（SSPs）下将物种分布预测到本世纪中叶（2041–2060）：SSP1-2.6、SSP2-4.5、SSP3-7.0 和 SSP5-8.5。为了避免与未来预测不可用或高度不确定的变量（如详细的土地利用或土壤属性）相关的外推风险，我们仅使用六个情景一致的变量重新训练所有模型：温度、降水、季节性、海拔和坡度。这些变量具有特征明确的未来预测，对淡水鱼类生理和栖息地结构至关重要，并且可以直接与 DAG 分析中识别的因果路径相关联。重新训练的模型确保预测是在训练数据的域内进行的，保持了稳健的性能（AUC ≥ 0.87）。

对于每个 SSP 和算法，我们生成了连续的适宜性地图，并总结了跨越两个正交不确定性维度的集合行为。结构不确定性定义为给定情景下四种建模算法预测的标准差，反映了算法假设和灵活性的差异（Lawler 等，2006；Thuiller 等，2019）。情景不确定性定义为给定算法下四种 SSP 路径预测的标准差，反映了排放和气候轨迹之间的分歧。我们使用方差划分来比较研究区域内这些不确定性来源的相对幅度，并识别保护建议对结构和情景不确定性均稳健的河流路段。

## 3. 结果

### 3.1 模型性能与整体拟合

四种核心算法在独立测试集上表现出稳健的预测性能。`Maxnet` 实现了最高的判别能力，AUC 为 0.909，敏感性为 0.896，表明其正确识别存在位置的能力很强。`随机森林` (RF) 在真实技巧统计量 (TSS = 0.699) 和特异性 (0.789) 方面表现出色，表明它在最小化假阳性方面特别有效。`GAM` 提供了相当的准确性 (AUC 0.897)，在复杂性和可解释性之间提供了平衡。`神经网络` (NN) 保持了可接受的判别力 (AUC 0.813)，但与其他方法相比，在预测背景样本时更为保守。跨算法的性能指标详细比较见表 2。这些模型的集合提供了栖息地适宜性的全面视图，利用了每种算法方法的优势，并突出了模型一致和不一致的区域。

**[插入表 2：Maxnet、RF、GAM 和神经网络的模型性能指标（AUC、TSS、敏感性、特异性、Boyce 指数）。]**

总体而言，判别指标超过了通常用于“良好”模型性能的阈值（AUC > 0.8），表明所选预测因子捕捉到了塑造鲫鱼分布的关键环境梯度（Elith 和 Leathwick，2009；Valavi 等，2021）。然而，与之前关于淡水入侵物种和 SDM 集合的研究一致（Rodríguez-Rey 等，2019；Hao 等，2020），我们将这些指标视为必要但不充分的条件。高 AUC 和 TSS 值可能与有偏的外推或错误指定的因果结构共存，特别是在受人类严重改造的河网中。因此，我们随后的因果推断和不确定性分解提供了必不可少的第二层评估，阐明了哪些预测因子可能是机制驱动因素，以及预测在何处仍然脆弱。

### 3.2 关键环境驱动因素与因果效应

多模型变量重要性分析表明，地形坡度和上游气候梯度是鲫鱼分布的一致驱动因素。`Maxnet` 和 `NN` 对平均海拔 (`dem_avg`) 和坡度变异性 (`slope_range`) 排名较高，而 `RF` 突出了河网拓扑指标，如汇流累积量 (`flow_acc`)。`GAM` 平滑项显示平均海拔对概率的非线性贡献最大。推断出的因果结构（如图 3 所示）解开了这些相关性。

**[插入图 3：共识因果 DAG（有向无环图），显示变量与物种出现之间的稳定依赖关系。]**

因果分析提供了比变量重要性排名更细致的视图。基于双重机器学习的平均处理效应 (ATE) 估算阐明了因果贡献，将其与单纯的关联区分开来（表 3）。城市建成区比例显示出对适宜性的显著正向因果效应 (ATE = 0.247 ± 0.026, p < 10⁻²⁰)，表明通常富含营养的低坡度城市水体对于这种耐受性物种来说仍然是有价值的栖息地。年降水量和上游粘土含量也是显著的正向因果驱动因素。重要的是，海拔范围和坡度变异性的正向效应证实了复杂地形中的河谷栖息地提供有利条件的假设，这可能是由于水文稳定性和多样化的微生境。

**[插入表 3：通过双重机器学习估算的关键环境驱动因素对物种出现的平均处理效应 (ATE)。]**

### 3.3 空间适宜性模式与不确定性

在当前气候条件下，高概率像素主要分布在长江中下游、江汉平原、珠江三角洲和松花江干流，如图 4 所示。这些区域对应于低海拔、高汇流累积量且水资源丰富的地区。预测标准差的空间均值为 0.099，表明模型一致性普遍较高。然而，方差分解揭示了一个显著的模式：结构不确定性（跨算法方差）平均约为情景不确定性（跨排放路径方差）的 4.7 倍。这表明，建模技术的选择比未来气候情景的选择给预测带来了更多的不确定性，凸显了集合建模的重要性。

**[插入图 4：当前栖息地适宜性（集合平均）地图及结构不确定性（跨模型标准差）的空间分布。]**

### 3.4 异质响应与局部敏感性

局部敏感性分析显示，物种对环境变量的响应存在显著的空间梯度。南部河流路段对坡度变异性表现出最高的平均敏感性，这意味着这些地区地形的微小变化可能对适宜性产生巨大影响。相比之下，中部和北部路段对季风降水变异性表现出更强的负向响应。CATE 分析（图 5）显示，虽然大多数河流路段对城市干预表现出较小的正向响应（与 ATE 一致），但在干旱的西北部出现了一个明显的负向响应尾部。这表明在缺水地区，城市化及相关的取水可能是有害的，这与在水资源丰富地区看到的正向效应形成对比。

**[插入图 5：关键驱动因素（如城市土地覆盖）的条件平均处理效应 (CATE) 空间图，突出正向与负向响应的区域。]**

### 3.5 未来情景下的分布变化

基于 CMIP6 SSP 情景对 2041–2060 年期间的预测显示，栖息地适宜性呈一致下降趋势。与当前分布相比，`Maxnet` 平均适宜性下降了约 25%，`GAM` 下降了约 50%。空间预测图表明，当前的热点河流路段将显著萎缩，北部冰雪补给的河流将受到热量增加和水文状况改变的不利影响。尽管在幅度上存在结构不确定性，但在不同的 SSP 和算法中，高适宜性路段的下降方向高度一致，为该物种的气候脆弱性提供了强有力的信号。图 6 直观展示了不同情景下的这些未来预测。

**[插入图 6：四种 SSP 情景（SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5）下的未来栖息地适宜性预测（2041–2060）。]**

## 4. 讨论

我们将因果发现与物种分布建模相结合，解决了相关性框架的一个根本局限：无法区分直接驱动因素与混杂关联（Pearl，2000；Peters 等，2017）。稳定的因果网络揭示了分层依赖关系，其中地形影响气候，进而塑造土地覆盖和土壤属性以构建栖息地适宜性。这与传统的变量重要性排名形成对比，后者通常混淆了直接和间接效应（Elith 和 Leathwick，2009；Molnar，2020）。例如，虽然海拔（`dem_avg`）在所有 SDM 中都显得非常重要，但因果分析表明，它主要通过调节季节性温度和降水起作用，而不是作为直接的生理限制。这种机制上的清晰度使得针对性干预成为可能；管理者可以解决特定的气候-地形相互作用，如维持冷空气汇聚或管理地形降水输入，而不是专注于“高海拔保护”。因果框架使我们能够剥离相关性层，识别真正控制生物多样性分布的杠杆。

河网一致性的实施改变了水生 SDM 的准确性和可解释性。大多数淡水 SDM 使用叠加在河流线上的陆地网格单元来近似暴露度，忽略了上游累积和连通性（Domisch 等，2015；Karger 等，2017）。通过采用 EarthEnv-Streams 和相关的水文环境图层，我们捕捉到了真实的水生暴露度：一个河流像素的环境整合了整个上游流域，而不仅仅是其局部的 1 km² 单元（Barbarossa 等，2018；Linke 等，2019）。这一区别至关重要——`flow_acc` 和 `flow_length` 在顶级预测因子中名列前茅，反映了物种对陆地网格所忽略的流域尺度干扰的敏感性。我们的结果可推广到任何河流生物，并可纳入描述流量状况改变、连通性丧失和生态网络结构的新兴网络指标（Tonkin 等，2018）。它们强调，对于河流物种而言，“环境”不是空间中的一个点，而是上游景观的累积遗产。

另一个惊人的发现是，跨模型方差（结构不确定性）使跨情景方差相形见绌，这表明在这种情况下，算法选择对预测不确定性的主导作用超过了排放路径。这与最近要求优先量化结构不确定性的呼吁一致，但与在没有模型集合的情况下报告情景范围的常见做法形成对比（Lawler 等，2006；Thuiller，2004；Thuiller 等，2019）。我们的结果表明，为了进行稳健的保护规划，管理者应采用多模型集合，而不是依赖单一的“最佳”算法，并专注于空间共识区域，同时将分歧区域视为需要在深度不确定性下进行适应性管理的区域。“情景可用变量再训练”策略在这方面被证明是必不可少的。通过确保当前和未来条件下的特征空间相同，我们避免了外推到未观测到的协变量空间——这是 SDM 转移到未来气候时普遍存在但经常被忽视的风险（Araújo 和 New，2007；Yates 等，2018）。

CATE 框架通过量化环境变化在*何处*产生最大影响，进一步通过因果推断实现可操作化。中海拔支流中的高 CATE 区域表明了恢复投资的“杠杆点”，而负 CATE 的低地城市流域可能需要替代策略，如缓解城市热岛效应或雨水管理。这种空间定位超越了忽略效应异质性的统一保护方法，并与关于淡水生物多样性空间明确气候适应规划的新兴工作产生共鸣（Reid 等，2019；Tonkin 等，2018）。然而，局限性仍然存在。因果发现算法假设无环性和因果充分性。虽然 bootstrap 稳定性过滤减轻了虚假边缘，但未测量变量（如生物相互作用或扩散障碍）残留的混杂因素仍然可能存在。未来的工作应结合先验知识，整合捕捉生物关联的联合物种分布模型（Tikhonov 等，2020），并测试替代算法以进一步完善这些因果图。

## 5. 结论

在本研究中，我们证明了将相关性 SDM 嵌入因果推断和河网框架中，可以对气候变化下的淡水鱼类分布产生更具可解释性和管理相关性的预测。通过将高分辨率水文环境图层与因果发现、用于平均处理效应的双重机器学习以及基于 CATE 的局部敏感性分析相结合，我们识别了一小组稳健的环境驱动因素和空间杠杆点，尽管存在巨大的结构不确定性，这些因素在不同算法和排放情景中仍保持一致。

我们的结果强调，对于河流物种，算法选择和网络结构的明确表示可以主导预测不确定性，这意味着保护规划应优先考虑多模型集合、情景可用预测因子以及对结构方差的透明核算。所提出的工作流程可直接转移到存在可比数据的其他淡水分类群和地区，并随着这些数据的可用性扩展到联合物种分布模型和动态连通性指标。通过从“预测并描述”转变为“推断并干预”范式，我们的方法为将大规模生物多样性预测与空间针对性适应行动联系起来提供了一个具体模板。

## 参考文献

1.  Araújo, M. B., & Guisan, A. (2006). Five (or so) challenges for species distribution modelling. *Journal of Biogeography*, 33(10), 1677–1688.
2.  Araújo, M. B., & New, M. (2007). Ensemble forecasting of species distributions. *Trends in Ecology & Evolution*, 22(1), 42–47.
3.  Athey, S., & Imbens, G. (2016). Recursive partitioning for heterogeneous causal effects. *Proceedings of the National Academy of Sciences*, 113(27), 7353–7360.
4.  Barbarossa, V., et al. (2018). FLO1K, global maps of mean, maximum and minimum annual streamflow at 1 km resolution from 1960 through 2015. *Scientific Data*, 5, 180052.
5.  Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32.
6.  Chernozhukov, V., et al. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68.
7.  Comte, L., & Olden, J. D. (2021). Evidence for dispersal syndromes in freshwater fishes. *Proceedings of the Royal Society B*, 288(1951), 20210223.
8.  Domisch, S., Amatulli, G., & Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution. *Scientific Data*, 2, 150073.
9.  Domisch, S., et al. (2019). Social equity shapes zone-selection: Balancing aquatic biodiversity representation and ecosystem services delivery. *Scientific Reports*, 9, 3082.
10. Dormann, C. F., et al. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, 36(1), 27–46.
11. Dudgeon, D., et al. (2006). Freshwater biodiversity: Importance, threats, status and conservation challenges. *Biological Reviews*, 81(2), 163–182.
12. Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics*, 40, 677–697.
13. Elith, J., et al. (2011). A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions*, 17(1), 43–57.
14. Filipe, A. F., et al. (2013). Biogeography of Iberian freshwater fishes revisited: The roles of historical versus contemporary constraints. *Journal of Biogeography*, 36(11), 2096–2110.
15. Guisan, A., & Zimmermann, N. E. (2000). Predictive habitat distribution models in ecology. *Ecological Modelling*, 135(2-3), 147–186.
16. Guisan, A., et al. (2013). Predicting species distributions for conservation decisions. *Ecology Letters*, 16(12), 1424–1435.
17. Hao, T., et al. (2020). Testing whether ensemble modelling is advantageous for maximising predictive performance of species distribution models. *Ecography*, 43(4), 549–558.
18. Hastie, T., & Tibshirani, R. (1986). Generalized additive models. *Statistical Science*, 1(3), 297–318.
19. He, Q., & Silliman, B. R. (2019). Climate change, human impacts, and coastal ecosystems in the Anthropocene. *Current Biology*, 29(19), R1021–R1035.
20. Hirzel, A. H., et al. (2002). Ecological-niche factor analysis: How to compute suitability maps without absence data? *Ecology*, 83(7), 2027–2036.
21. Kalisch, M., & Bühlmann, P. (2007). Estimating high-dimensional directed acyclic graphs with the PC-algorithm. *Journal of Machine Learning Research*, 8, 613–658.
22. Karger, D. N., et al. (2017). Climatologies at high resolution for the earth's land surface areas. *Scientific Data*, 4, 170122.
23. Lawler, J. J., et al. (2006). Predicting climate-induced range shifts: Model differences and model reliability. *Global Change Biology*, 12(8), 1568–1584.
24. Linke, S., et al. (2019). Global hydro-environmental sub-basin and river reach characteristics at high spatial resolution. *Scientific Data*, 6, 283.
25. Ma, B., et al. (2020). Machine learning for causal inference in biology. *Nature Methods*, 17, 1129–1130.
26. Merow, C., et al. (2013). A practical guide to MaxEnt for modeling species' distributions: What it does, and why inputs and settings matter. *Ecography*, 36(10), 1058–1069.
27. Molnar, C. (2020). *Interpretable machine learning: A guide for making black box models explainable*. Lulu.com.
28. Olden, J. D., et al. (2010). Conservation biogeography of freshwater fishes: Recent progress and future challenges. *Diversity and Distributions*, 16(3), 496–513.
29. Pearl, J. (2000). *Causality: Models, Reasoning, and Inference*. Cambridge University Press.
30. Pearl, J. (2009). Causal inference in statistics: An overview. *Statistics Surveys*, 3, 96–146.
31. Pearson, R. G., & Dawson, T. P. (2003). Predicting the impacts of climate change on the geographical distribution of species: Criticisms, uncertainties, and case studies. *Global Ecology and Biogeography*, 12(5), 361–371.
32. Peters, J., et al. (2017). *Elements of Causal Inference: Foundations and Learning Algorithms*. MIT Press.
33. Peterson, A. T., et al. (2011). *Ecological Niches and Geographic Distributions*. Princeton University Press.
34. Phillips, S. J., et al. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3-4), 231–259.
35. Phillips, S. J., et al. (2017). Opening the black box: An open-source release of Maxent. *Ecography*, 40(7), 887–893.
36. Poff, N. L., et al. (1997). The natural flow regime. *BioScience*, 47(11), 769–784.
37. Qiao, H., et al. (2015). Theoretical underpinnings of ecological niche modeling: A conceptual review. *Ecological Modelling*, 312, 56–68.
38. Reid, A. J., et al. (2019). Emerging threats and persistent conservation challenges for freshwater biodiversity. *Biological Reviews*, 94(3), 849–873.
39. Schapire, R. E. (2003). The boosting approach to machine learning: An overview. *Nonlinear Estimation and Classification*, 171, 149–171.
40. Shapley, L. S. (1953). A value for n-person games. *Contributions to the Theory of Games*, 2, 307–317.
41. Spirtes, P., et al. (2000). *Causation, Prediction, and Search*. MIT Press.
42. Stockwell, D., & Peters, D. (1999). The GARP modelling system: Problems and solutions to automated spatial prediction. *International Journal of Geographical Information Science*, 13(2), 143–158.
43. Su, Y., et al. (2021). An overview of causal inference in ecology. *Ecological Indicators*, 129, 107960.
44. Thuiller, W. (2004). Patterns and uncertainties of species' range shifts under climate change. *Global Change Biology*, 10(12), 2020–2027.
45. Thuiller, W., et al. (2009). BIOMOD – a platform for ensemble forecasting of species distributions. *Ecography*, 32(3), 369–373.
46. Thuiller, W., et al. (2019). Uncertainty in ensembles of global biodiversity scenarios. *Nature Communications*, 10, 1446.
47. Tikhonov, G., et al. (2020). Joint species distribution modelling with the r-package Hmsc. *Methods in Ecology and Evolution*, 11(3), 442–447.
48. Tonkin, J. D., et al. (2018). Flow regime alteration degrades ecological networks in riparian ecosystems. *Nature Ecology & Evolution*, 2, 86–93.
49. Valavi, R., et al. (2021). predictive performance of presence-only species distribution models: a benchmark study with reproducible code. *Ecological Monographs*, 92(1), e01486.
50. Van der Laan, M. J., & Rose, S. (2011). *Targeted Learning: Causal Inference for Observational and Experimental Data*. Springer.
51. Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242.
52. Wisz, M. S., et al. (2008). Effects of sample size on the performance of species distribution models. *Diversity and Distributions*, 14(5), 763–773.
53. Rodríguez-Rey, M., Consuegra, S., Börger, L., & Garcia de Leaniz, C. (2019). Improving Species Distribution Modelling of freshwater invasive species for management applications. *PLoS ONE*, 14(6), e0217896.
54. WWF. (2020). *Living Planet Report 2020 – Bending the curve of biodiversity loss*. WWF.
55. Yates, K. L., et al. (2018). Outstanding challenges in the transferability of ecological models. *Trends in Ecology & Evolution*, 33(10), 790–802.