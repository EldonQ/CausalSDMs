# CAST: Causal Analysis and Screening Toolkit for Species Distribution Modeling

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Research_Article-orange.svg)]()

> **[中文文档 (Chinese Version)](README_CN.md)**

## 📖 Overview

**CAST** (Causal Analysis and Screening Toolkit) is a methodological framework designed to advance species distribution modeling (SDM) from **correlative prediction** to **mechanistic inference**. By integrating Causal Discovery (Bayesian Networks) with Double Machine Learning (DML), CAST explicitly separates true causal drivers from spurious correlations. 

This repository contains the official implementation of the CAST framework, including the complete analytical pipeline, from environmental data processing to the construction of Causal-Informed Multi-Layer Perceptrons (CI-MLP).

## ✨ Key Capabilities

- **Automated Causal Discovery**: Learns causal Directed Acyclic Graphs (DAGs) from observational data to identify hierarchical relationships among environmental drivers. 
- **Causal Effect Estimation**: Quantifies Average Treatment Effects (ATE) using Double Machine Learning, decoupling true ecological mechanisms from mere predictive importance.
- **Adaptive Feature Screening**: Provides a scoring algorithm to filter out redundant and spurious variables based on local DAG topology, causal significance, and predictive contribution.
- **Causal-Informed Neural Networks (CI-MLP)**: Embeds continuous structural knowledge (ATE-derived feature weights and DAG-guided interaction terms) into a deep learning architecture, significantly improving predictive generalization and preventing overfitting.
- **Spatially Explicit Causal Mapping**: Evaluates Conditional Average Treatment Effects (CATE) via Causal Forests to map heterogeneous environmental responses across geographical spaces.

## 📂 Repository Structure

```text
E:\CausalSDMs\
├── data-main/               # Environmental datasets (Topography, Climate, Hydrology)
├── scripts/EcoISEA3H/       # Core CAST analytical pipeline
│   ├── 01_data_preparation...  # Data aggregation and background sampling
│   ├── 02_main_cast_pipeline/  # DAG learning, ATE estimation & CAST screening
│   ├── 03_run_Eco_multi_sp/    # Batch CI-MLP modeling and baseline comparisons
│   └── plot/                   # Figure generation scripts (Nature-standard)
├── outputs/                 # Intermediate results, continuous weights, & tables
├── figures/                 # Final manuscript graphics
└── manuscript_mee_cn.md     # Manuscript drafts and discussions
```

## 🚀 Usage Guide

The CAST workflow operates sequentially. Key implementation scripts include:

1.  **Causal Graph Learning & ATE Estimation**
    Executes constraints-based structure learning and DML algorithms.
    *Script*: `scripts/EcoISEA3H/02_cast_screening_pipeline.R`
2.  **Adaptive Variable Screening**
    Computes component scores (DAG, ATE, RF) and extracts the final causal variable subset.
    *Script*: `scripts/EcoISEA3H/02_cast_screening_pipeline.R`
3.  **CI-MLP Network Training**
    Trains the Causal-Informed MLP using weighted inputs and DAG topology interactions, evaluating performance against standard machine learning models (FlatNN, RF, Maxent, BRT).
    *Script*: `scripts/EcoISEA3H/03_run_Eco_multi_species.R`
4.  **CATE Spatial Mapping**
    Generates high-resolution spatial models of Conditional Average Treatment Effects.
    *Script*: `scripts/EcoISEA3H/04_spatial_cate.R`

## 🛠️ Requirements & Dependencies

*   **R** (version >= 4.1.0)
*   **Causal Inference**: `bnlearn`, `pcalg`, `DoubleML`, `grf`
*   **Machine Learning**: `keras`, `tensorflow`, `randomForest`, `maxnet`, `gbm`
*   **Geospatial**: `terra`, `sf`
*   **Visualization**: `ggplot2`, `patchwork`, `igraph`, `ggraph`

## 📜 Citation

If you use CAST in your research, please cite our associated methodological paper:

> *Mechanism-driven species distribution modelling via causal structural learning and effect estimation.* (Under Review)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
