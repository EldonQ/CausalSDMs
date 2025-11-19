# 3. Results

本研究聚焦于构建因果驱动的河流网络栖息地模型，以回答大陆尺度上中华鲫（*Carassius auratus*）在多源环境梯度下的可适生空间分布及其驱动机制。以下结果严格沿研究问题展开，直接回应所采用的方法设置，为后续讨论提供定量依据。

## 3.1 模型性能与整体拟合

四套核心算法在 336 个存在记录与 1680 个背景点上表现稳定，`Maxnet` 模型获得最高 AUC（0.909）与灵敏度（0.896），`RF` 模型在 TSS（0.699）与特异度（0.789）上占优，`GAM` 模型提供与树模型一致的综合准确度（AUC 0.897），而 `NN` 虽保持可接受的判别力（AUC 0.813），但对背景样本更为保守（阈值≈4.3×10⁻⁴），提示直连层数尚未完全捕获非线性关系。[Table placeholder: output/08_model_evaluation/evaluation_summary.csv] 系统汇总了全部评估指标，原始 CSV 供后续绘图与复核。

四模型的 ROC 曲线与阈值响应在 `[Figure placeholder: figures/08_model_evaluation/performance_comparison.png]` 与 `[Figure placeholder: figures/08_model_evaluation/roc_curves.png]` 中展示，曲线面积与纵向误差条均采用 95% 置信区间表达，确保读者可直观比较模型差异。

## 3.2 关键环境驱动与因果效应

多模型变量重要度分析表明，地形坡度与上游气候梯度是跨算法的一致驱动因素：`Maxnet` 与 `NN` 将 `dem_avg`、`slope_range` 及最湿季温度（`hydro_wavg_08`）列入前三；`RF` 则突出河网拓扑（`flow_acc`、`flow_length`）与最湿季降水（`hydro_wavg_18`）；`GAM` 的平滑项显示平均高程对概率贡献最大（标准化 SHAP=0.206），其次是季节性温度（`hydro_wavg_08`）与最湿月降水（`hydro_wavg_13`）。[Table placeholder: output/09_variable_importance/importance_summary.csv] 与 `[Figure placeholder: figures/09_variable_importance/importance_heatmap.png]`、`[Figure placeholder: figures/09_variable_importance/importance_lollipop_by_model.png]` 记录了模型间的一致性与差异性。

基于 Double Machine Learning 的平均处理效应（ATE）估计进一步厘清了土地利用与水文条件的因果贡献：城市建成区占比（`lc_wavg_09`，ATE=0.247±0.026，p<10⁻²⁰）显著提高适宜性，说明低坡缓流的城市水体对中华鲫仍具栖息价值；年降水（`hydro_wavg_12`，ATE=0.146±0.009）与上游黏土含量（`soil_wavg_05`，ATE=0.143±0.014）同样正向显著，而高程范围（`dem_range`）与坡度变幅（`slope_range`）的正效应印证了在复杂地形中存在高适宜性谷地的假设。[Table placeholder: output/14_causal/ate_all_variables.csv] 保存了全部变量的估计值。[Figure placeholder: figures/14_causal/ate_all_variables_forest.png] 与 `[Figure placeholder: figures/14_causal/dag_core_pathways_accurate.png]` 展示了显著边的影响路径与结构图。

## 3.3 空间适宜性格局与不确定性

基准情景下，四模型在 1.86 百万个河网像元上的预测均值介于 0.099（`NN`）至 0.201（`RF`）之间，中位数低于 0.10 而 90 分位数攀升至 0.49–0.67，指向“低适宜背景 + 高适宜河段斑块”的格局。[Figure placeholder: figures/11_prediction_maps/prediction_rf.png]、`[Figure placeholder: figures/11_prediction_maps/prediction_maxnet.png]`、`[Figure placeholder: figures/11_prediction_maps/prediction_gam.png]` 与 `[Figure placeholder: figures/11_prediction_maps/prediction_nn.png]` 显示，高概率栅格主要沿长江中下游、珠江三角洲与松花江干流分布，两湖平原及太湖流域呈连续热点，与历史记录集中区吻合。

模型间一致性通过标准差和一致性指数量化：预测标准差空间均值为 0.099（90 分位 0.264），模型一致性均值 0.780，90 分位达 0.995，说明高概率河段的跨模型共识度极高；结合未来情景的方差分解结果可以发现，结构性不确定性（跨模型标准差）在全国范围内显著高于情景不确定性（跨 SSP 标准差），平均高约 4.7 倍，而低概率地区的不确定性总体较低（`[Figure placeholder: figures/12_uncertainty/uncertainty_map.png]`、`[Figure placeholder: figures/12_uncertainty/model_agreement.png]`）。[Table placeholder: output/12_uncertainty/uncertainty_summary.csv] 记录了不确定性分布，用于后续统计检验。

此外，基于因果森林估计的个体化处理效应（CATE）揭示局地异质性：`cate_summary.csv`（均值 0.013，中位数 0.013，90 分位 0.064）显示大部分河段对建成区干预保持小幅正响应，西北干旱区出现负响应尾部（最小值 -0.146），提示未来环境管理应关注城市扩张对干流热点的影响。[Figure placeholder: figures/11_prediction_maps/cate_map.png] 提供空间化分布。

## 3.4 异质响应与局地敏感度

局地敏感度分析按纬向分区呈现显著梯度。[Table placeholder: output/17_local_sensitivity/sensitivity_summary.csv] 表明南方河段对坡度变幅的平均敏感度最高（0.022±0.019），中北部则对季风降水（`hydro_wavg_08`）表现出更强的负响应（北方平均 -0.037）。耕地占比（`lc_wavg_07`）在全国范围内均呈轻微正向敏感度，而土壤黏土（`soil_wavg_05`）与砂含量（`soil_wavg_03`）在南北方表现为负敏感度，映射出沉积物负荷与底质对栖息地质量的约束。[Figure placeholder: figures/17_local_sensitivity/sensitivity_heatmap.png] 与 `[Figure placeholder: figures/17_local_sensitivity/sensitivity_violin.png]` 为深入解析提供视觉支持。

## 3.5 未来情景下的分布变化

在仅使用 6 个跨时段可用的 Bioclim 与地形变量（`bio01`、`bio12`、`bio04`、`bio15`、`dem_avg`、`slope_avg`）重新训练四套模型并在独立测试集上保持较高判别力（AUC≥0.87）的前提下，基于 CMIP6 SSP 情景的气候投影（截至 2060 年代）显示四模型的一致衰减趋势：相较当前分布，`Maxnet` 的平均适宜性在 SSP126/245/370/585 下分别下降 25.4%、24.9%、24.8% 与 24.1%，`GAM` 均值下降约 50%（0.158→0.078），`RF` 与 `NN` 的均值分别降至 0.029 与 0.045，说明未来河网整体适宜性收缩但仍保留局地高值核心。[Table placeholder: output/15_future_env/prediction_trends_all_models.csv] 与 `[Figure placeholder: figures/15_future_env/habitat_trends_all_models.png]` 汇总模型间差异，而各情景下的空间预测图（如 `[Figure placeholder: figures/15_future_env/SSP585/prediction_rf.png]`）显示热点河段缩小、北方冰雪补给河段受热量增加影响更显著。

为量化驱动背景，[Table placeholder: output/15_future_env/future_bioc_statistics.csv] 提供 Bioclim 变量的情景统计，揭示 SSP585 下年均温升高 1.0–1.5°C、极端降水加剧，与适宜性下降的幅度吻合；同时，未来评估结果表明，高适宜性河段在不同 SSP 与算法下的下降方向高度一致，而情景间差异相对较小。未来评估文件夹保留了全部预测栅格（1200 dpi PNG 与 SVG）及 `.tif` 数据，便于后续高分辨率制图与不确定性扩展分析。
