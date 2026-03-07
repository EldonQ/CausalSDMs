# CAST: 物种分布模型因果分析与变量筛选工具包

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Research_Article-orange.svg)]()

> **[English Version (英文原版)](README.md)**

## 📖 简介

**CAST (Causal Analysis and Screening Toolkit)** 是一个新的方法学框架，旨在将物种分布模型（SDMs）从**相关性预测**推进至**因果机制推断**。CAST 通过整合因果发现算法（有向无环图）与双重机器学习（DML），显式地将真实的生态因果驱动因子与虚假的统计相关性进行解耦。

本仓库包含了 CAST 框架的官方实现代码，涵盖从环境数据处理、因果推断到构建“因果启发的神经网络（CI-MLP）”的完整分析流水线。

## ✨ 核心功能

- **自动化因果发现**：基于观测数据学习环境因子的因果有向无环图 (DAG)，提炼环境变量间的层级与依赖关系。
- **因果效应估计**：利用双重机器学习（Double Machine Learning）量化环境变量的平均处理效应 (ATE)，将真实的生态机制与纯粹的预测重要性剥离。
- **自适应特征筛选**：提供一套综合评分算法，基于局部 DAG 拓扑结构、因果效应显著性及预测重要性，剔除冗余和伪相关变量。
- **因果启发的神经网络 (CI-MLP)**：将连续的因果结构先验（基于 ATE 算得的特征权重、基于 DAG 提取的交互项）嵌入深度学习架构中，显著提升模型的泛化能力并有效防止过拟合。
- **解译空间异质性**：通过因果森林 (Causal Forests) 评估条件平均处理效应 (CATE)，绘制环境变量驱动力的空间异质性地图。

## 📂 目录结构

```text
E:\CausalSDMs\
├── data-main/               # 环境基础数据 (地形、气候、水文数据集)
├── scripts/EcoISEA3H/       # CAST 核心分析代码库
│   ├── 01_data_preparation...  # 数据聚合与背景点采样
│   ├── 02_main_cast_pipeline/  # DAG 学习、ATE 估计与 CAST 变量筛选
│   ├── 03_run_Eco_multi_sp/    # 模型训练 (CI-MLP 构建) 及基线模型对比
│   └── plot/                   # 结果图表渲染代码 (对标 Nature 格式)
├── outputs/                 # 中间结果、连续权重系数及表格计算结果
├── figures/                 # 论文正式发布的最终图表
└── manuscript_mee_cn.md     # 论文中文底稿及修订历史
```

## 🚀 核心工作流

CAST 提供结构化的模块代码，核心执行步骤如下：

1.  **因果图学习与 ATE 估计**
    运行基于约束和评分的拓扑学习，并执行 DML 算法以估算因果效应。
    *脚本*: `scripts/EcoISEA3H/02_cast_screening_pipeline.R`
2.  **自适应因果变量筛选**
    计算 DAG、ATE 及 RF 三大组件得分，输出用于最终生态建模的纯净因果变量子集。
    *脚本*: `scripts/EcoISEA3H/02_cast_screening_pipeline.R`
3.  **CI-MLP 神经网络训练**
    构建因果启发的感知机网络（导入 ATE 权重与 DAG 交互项），并与 FlatNN、RF、Maxent 等基线模型执行标准化对比评估。
    *脚本*: `scripts/EcoISEA3H/03_run_Eco_multi_species.R`
4.  **因果效应空间制图 (CATE)**
    构建高分辨率的条件平均处理效应空间分布模型，辅助局地保护决策。
    *脚本*: `scripts/EcoISEA3H/04_spatial_cate.R`

## 🛠️ 环境与依赖

*   **R 语言环境** (版本 >= 4.1.0)
*   **因果推断核心库**: `bnlearn`, `pcalg`, `DoubleML`, `grf`
*   **机器学习算法库**: `keras`, `tensorflow`, `randomForest`, `maxnet`, `gbm`
*   **空间分析**: `terra`, `sf`
*   **图表可视化**: `ggplot2`, `patchwork`, `igraph`, `ggraph`

## 📜 引用

如果您在学术研究中使用了 CAST 分析框架，请引用我们关联的方法学论文：

> *Mechanism-driven species distribution modelling via causal structural learning and effect estimation.* (Under Review)

## 📄 许可证

本项目采用 MIT 许可证开源 - 详情请参阅 [LICENSE](LICENSE) 文件。
