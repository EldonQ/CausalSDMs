## 项目概述

本项目使用四种机器学习算法(Maxnet, Random Forest, GAM, Neural Network)对Carassius auratus（鲫鱼）在中国淡水系统中的分布进行建模和预测,旨在发表于Nature级别期刊。

**研究物种**: Carassius auratus (鲫鱼/Crucian Carp)
**研究区域**: 中国淡水系统(水网)  
**环境数据**: EarthEnv Streams (324层淡水环境变量，以选取有意义的173个环境变量,1km分辨率)

---

## 脚本说明

所有脚本按照执行顺序编号,高度模块化,可独立运行。每个脚本都会生成详细的日志文件。

### 数据准备阶段

#### 01_data_preparation_NEW.R
- **功能**: 物种出现数据的清洗和质量控制
- **输入**: `Carassius_auratusOCC/Carassius_auratus.csv` (GBIF数据)
- **输出**: 
  - `output/01_data_preparation/species_occurrence_cleaned.csv` - 清洗后的物种数据
  - `output/01_data_preparation/species_occurrence_cleaned.shp` - Shapefile格式
- **主要步骤**:
  1. 读取GBIF数据
  2. 标准化数据格式
  3. 限制到中国境内(使用中国边界)
  4. 去除重复和异常坐标
  5. 空间稀疏化(约10km网格)

#### 02_env_extraction.R
- **功能**: 从netCDF环境变量文件中提取物种出现点的环境数据
- **输入**: 
  - `output/01_species_occurrence_cleaned.csv`
  - `earthenvstreams/*.nc`
- **输出**: 
  - `output/02_occurrence_with_env.csv` - 带环境变量的物种数据
  - `output/02_extracted_variables.csv` - 提取的变量列表
- **主要步骤**:
  1. 读取324层环境变量
  2. 提取关键变量(水文气候、地形、流量、土地覆盖、土壤)
  3. 应用单位转换
  4. 质量检查和缺失值处理

#### 03_background_points.R
- **功能**: 生成背景点(伪不存在点)
- **输入**: `output/02_occurrence_with_env.csv`
- **输出**: 
  - `output/03_background_points.csv` - 背景点数据
  - `output/03_combined_presence_absence.csv` - 合并的出现点和背景点
- **主要步骤**:
  1. 使用淡水掩膜限定背景点范围
  2. 生成10倍于出现点数量的背景点
  3. 提取背景点的环境变量
  4. 合并出现点和背景点数据

### 数据分析阶段

#### 04_collinearity_analysis.R
- **功能**: 环境变量相关性分析和多重共线性诊断
- **输入**: `output/03_combined_presence_absence.csv`
- **输出**: 
  - `output/04_correlation_matrix.csv` - 相关系数矩阵
  - `output/04_selected_variables.csv` - 筛选后的变量列表
  - `output/04_collinearity_removed.csv` - 最终建模数据集
  - `figures/04_correlation_heatmap.png` - 相关性热图
- **主要步骤**:
  1. 计算Pearson相关系数矩阵
  2. VIF(方差膨胀因子)分析
  3. 基于相关系数的逐步筛选(|r| < 0.8)
  4. 生成诊断图表

### 模型训练阶段

#### 05_model_maxnet.R
- **功能**: Maxnet(Maximum Entropy)模型训练和评估
- **输入**: `output/04_collinearity_removed.csv`
- **输出**: 
  - `output/05_maxnet_model.rds` - 训练好的模型
  - `output/05_maxnet_predictions.csv` - 预测结果
  - `output/05_maxnet_evaluation.csv` - 性能评估(AUC, TSS, Kappa)
  - `output/05_maxnet_variable_importance.csv` - 变量重要性
- **模型参数**: 
  - regmult = 1.0
  - feature classes = lqph (linear, quadratic, product, hinge)

#### 05b_model_nn.R
- **功能**: 神经网络(Neural Network)模型训练和评估
- **输入**: `output/04_collinearity_removed.csv`
- **输出**: 
  - 训练好的神经网络模型
  - 预测结果和性能评估
  - 变量重要性

#### 06_model_rf.R
- **功能**: 随机森林(Random Forest)模型训练和评估
- **输入**: `output/04_collinearity_removed.csv`
- **输出**: 
  - `output/06_rf_model.rds` - 训练好的模型
  - `output/06_rf_predictions.csv` - 预测结果
  - `output/06_rf_evaluation.csv` - 性能评估
  - `output/06_rf_variable_importance.csv` - 变量重要性
- **模型参数**: 
  - ntree = 500
  - mtry = sqrt(n_variables)
  - 平衡采样处理类别不平衡

#### 07_model_gam.R
- **功能**: 广义加性模型(GAM)训练和评估
- **输入**: `output/04_collinearity_removed.csv`
- **输出**: 
  - `output/07_gam_model.rds` - 训练好的模型
  - `output/07_gam_predictions.csv` - 预测结果
  - `output/07_gam_evaluation.csv` - 性能评估
  - `output/07_gam_variable_importance.csv` - 变量重要性
- **模型参数**: 
  - family = binomial(logit)
  - 平滑项 k = 5
  - method = fREML (fast REML)

#### 08_model_evaluation.R
- **功能**: 比较多个模型的性能并生成集成预测
- **输入**: 所有模型的评估和预测结果
- **输出**: 
  - `output/08_model_evaluation/model_comparison.csv` - 模型性能对比表
  - `output/08_model_evaluation/ensemble_predictions.csv` - 集成预测结果
  - `figures/08_model_evaluation/model_comparison_barplot.png` - 性能对比图
  - `figures/08_model_evaluation/roc_curves_comparison.png` - ROC曲线对比
- **集成方法**: 
  - 简单平均
  - 基于AUC的加权平均

### 结果可视化阶段

#### 09_variable_importance_viz.R
- **功能**: 变量重要性可视化(小提琴图、热图)
- **输入**: 三个模型的变量重要性结果
- **输出**: 
  - `figures/09_variable_importance_violin.png` - 小提琴图
  - `figures/09_variable_importance_top20.png` - Top 20变量条形图
  - `figures/09_variable_importance_heatmap.png` - 热图
- **分析内容**: 
  - 标准化变量重要性
  - 模型间一致性分析
  - Top变量识别

#### 10_response_curves.R
- **功能**: 生成环境变量响应曲线和偏依赖图
- **输入**: 训练好的模型和建模数据
- **输出**: 
  - `figures/10_response_curves_*.png` - 各模型响应曲线
  - `figures/10_response_curves_comparison.png` - 模型对比响应曲线
  - `output/10_response_data.csv` - 响应曲线原始数据
- **分析内容**: 
  - Top 10重要变量的响应曲线
  - 单变量边际效应
  - 生态阈值识别

#### 11_current_prediction_maps.R
- **功能**: 生成当前气候条件下的空间分布预测地图
- **输入**: 训练好的模型和环境变量栅格
- **输出**: 
  - `output/11_prediction_map_*.tif` - 预测栅格(GeoTIFF)
  - `figures/11_prediction_map_*.png` - 可视化地图
- **注意事项**: 
  - 使用降采样策略以节省内存
  - 生成Maxnet、RF、GAM和集成四种预测地图

#### 12_uncertainty_map.R
- **功能**: 生成模型预测不确定性地图
- **输入**: 三个模型的预测栅格
- **输出**: 
  - `output/12_uncertainty_map_sd.tif` - 标准差地图
  - `output/12_uncertainty_map_cv.tif` - 变异系数地图
  - `output/12_model_agreement_map.tif` - 模型一致性地图
  - `figures/12_*.png` - 可视化地图
- **不确定性指标**: 
  - 标准差(SD) - 绝对差异
  - 变异系数(CV) - 相对差异
  - 模型一致性 - 预测存在的模型数量

#### 13_study_area_maps.R
- **功能**: 生成研究区域和物种分布概况图
- **输入**: 清洗后的物种数据和背景点
- **输出**: 
  - `figures/13_study_area_map.png` - 研究区域总览
  - `figures/13_species_occurrence_map.png` - 物种出现点地图
  - `figures/13_combined_presence_absence_map.png` - 出现点+背景点地图
- **可视化内容**: 
  - 地理范围界定
  - 数据来源区分
  - 采样偏差评估

---

## 执行顺序

**严格按照脚本编号顺序执行**,每个脚本的输出是下一个脚本的输入:

```
01_data_preparation_NEW.R
  ↓
02_env_extraction.R
  ↓
03_background_points_NEW.R
  ↓
04_collinearity_analysis.R
  ↓
├─ 05_model_maxnet.R
├─ 05b_model_nn.R
├─ 06_model_rf.R
└─ 07_model_gam.R
  ↓
08_model_evaluation.R
  ↓
├─ 09_variable_importance_viz.R
├─ 10_response_curves.R
├─ 11_current_prediction_maps.R
├─ 12_uncertainty_map.R
├─ 13_study_area_maps.R
├─ 14_background_points_map.R
├─ 15_future_env_projection.R
├─ 16_future_prediction.R
└─ 17_local_sensitivity_analysis.R
```

---

## 运行方法

### 方法1: 单独运行每个脚本
```R
# 在R或RStudio中
source("scripts/01_data_preparation_NEW.R")
source("scripts/02_env_extraction.R")
# ... 依次执行
```

### 方法2: 批量运行
```bash
# 在命令行中
cd E:/SDM01
Rscript scripts/01_data_preparation_NEW.R
Rscript scripts/02_env_extraction.R
# ... 依次执行
```

---

## 输出文件结构

```
E:/SDM01/
├── output/              # 所有中间和最终数据结果
│   ├── 01_*.csv        # 数据准备阶段输出
│   ├── 02_*.csv        # 环境提取输出
│   ├── 03_*.csv        # 背景点输出
│   ├── 04_*.csv        # 共线性分析输出
│   ├── 05_*.rds/csv    # Maxnet模型输出
│   ├── 06_*.rds/csv    # RF模型输出
│   ├── 07_*.rds/csv    # GAM模型输出
│   ├── 08_*.csv        # 模型比较输出
│   ├── 09-13_*.csv     # 可视化相关数据
│   └── *_log.txt       # 各脚本的日志文件
│
└── figures/             # 所有科研图表(1200 dpi PNG)
    ├── 04_*.png        # 相关性热图
    ├── 08_*.png        # 模型比较图
    ├── 09_*.png        # 变量重要性图
    ├── 10_*.png        # 响应曲线
    ├── 11_*.png        # 预测分布图
    ├── 12_*.png        # 不确定性图
    └── 13_*.png        # 研究区域图
```

---

## 系统要求

### R版本
- R >= 4.0.0

### 必需的R包
```R
# 数据处理
tidyverse, data.table

# 空间分析
sf, sp, raster, ncdf4, rgdal

# 地图可视化
rnaturalearth, rnaturalearthdata, viridis, ggplot2, ggspatial

# 坐标清洗
CoordinateCleaner

# 模型包
maxnet, randomForest, mgcv, dismo, caret, pROC, nnet, neuralnet

# 其他
gridExtra, corrplot, car, usdm, pdp
```

### 硬件要求
- **内存**: 至少16GB RAM(推荐32GB)
- **存储**: 至少50GB可用空间
- **处理器**: 多核CPU(模型训练会自动使用并行计算)

---

## 注意事项

1. **数据完整性**: 确保`earthenvstreams/`文件夹包含所有必需的netCDF文件
2. **内存管理**: 脚本11会使用降采样策略,如果内存充足可以调整`downsample_factor`参数
3. **计算时间**: 
   - 脚本02(环境提取): 10-30分钟
   - 脚本05-07(模型训练): 每个5-20分钟
   - 脚本11(空间预测): 30-60分钟
4. **数据保留**: 所有中间数据都会保存,便于后续修改和重新绘图
5. **图表质量**: 所有图表均为1200 dpi,符合Nature期刊要求

---

## 引用数据源

### 环境数据
Domisch, S., Amatulli, G., Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1km resolution. Scientific Data, 2, 150073.

### 相关方法
- Maxnet: Phillips et al. (2017)
- Random Forest: Breiman (2001)
- GAM: Wood (2017)

---

## 联系方式

如有问题,请查看各脚本生成的日志文件(`output/*_log.txt`)进行故障排除。

---

**创建日期**: 2025-10-02  
**更新日期**: 2025-10-19  
**项目目标**: Nature级别期刊发表  
**研究主题**: 中国淡水系统中鲫鱼(Carassius auratus)的物种分布模型

> d <- read.csv('E:/SDM01/output/04_collinearity/selected_variables.csv'); lib$

变量类型统计:
                  type count
1 水文气候 (hydro_avg)     2
2      温度 (tmax_avg)     1
3      降水 (prec_sum)     1
4    土地覆盖 (lc_avg)    11
5      土壤 (soil_avg)     7
6      地质 (geo_wsum)    32

总计: 54 个变量