# CausalSDMs 🌊🐟

### Mechanism-Driven Freshwater Fish Distribution Modeling via Causal Inference

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Research_Prototype-orange.svg)]()

> **[中文说明 (Chinese Version)](README_CN.md)**

## 📖 Overview

**CausalSDMs** is a novel research framework designed to advance freshwater biodiversity forecasting from **correlative prediction** to **mechanistic inference**. By integrating **Causal Discovery** (Bayesian Networks) with traditional **Species Distribution Models (SDMs)**, this project addresses the challenge of distinguishing true causal drivers from spurious associations in complex river networks.

Applied to freshwater fish distributions across China, this framework leverages **EarthEnv-Streams** (upstream-weighted environmental layers) to capture the true aquatic exposure of riverine organisms.

## ✨ Key Features

- **🌊 River Network Consistency**: Utilizes 1-km resolution upstream-weighted variables (flow accumulation, hydro-climatic aggregation) rather than terrestrial grid overlays.
- **🔍 Causal Structure Learning**: Implements constraint-based (**PC**) and score-based (**Hill-Climbing**) algorithms to discover stable environmental dependency networks (DAGs).
- **📉 Causal Variable Selection**: Reduces predictor dimensionality by ~38% while maintaining or improving model performance (AUC/TSS) by filtering out non-causal confounders.
- **🎯 Heterogeneous Treatment Effects**: Estimates Conditional Average Treatment Effects (**CATE**) via Causal Forests to identify spatially explicit high-leverage zones for conservation.
- **🔮 Robust Future Projections**: distinguishing **Structural Uncertainty** (algorithmic) from **Scenario Uncertainty** (emission pathways).

## 📂 Repository Structure

```text
E:\CausalSDMs\
├── data-main/               # Core environmental datasets (Vector/Raster)
├── scripts/                 # Analysis pipeline (numbered by workflow step)
│   ├── 01-04_...           # Data preparation & collinearity analysis
│   ├── 05-08_...           # SDM training (Maxent, RF, GAM, NN) & Evaluation
│   ├── 09-13_...           # Variable importance, Response curves, Maps
│   ├── 14_causal...        # Causal discovery & ATE estimation
│   └── 15_future...        # Future climate projections & retraining
├── figures/                 # Generated plots (Nature-standard formatting)
├── output/                  # Model objects, intermediate CSVs, and Rasters
├── manuscript_*.md          # Drafts and core content for publication
└── README.md                # Project documentation
```

## 🚀 Workflow

The analysis pipeline is organized into numbered scripts for reproducibility:

1.  **Data Preparation**: Cleaning occurrence records and extracting EarthEnv-Streams variables (`01_data_preparation_NEW.R`).
2.  **Base Modeling**: Training four algorithms (Maxent, Random Forest, GAM, Neural Network) on the full variable set (`05_model_maxnet.R` to `07_model_gam.R`).
3.  **Causal Discovery**:
    *   Infer Causal DAGs using PC and Hill-Climbing algorithms (`14_causal_discovery.R`).
    *   Estimate Average Treatment Effects (ATE) via Double Machine Learning (`14c_batch_ate_estimation.R`).
4.  **Causal-Informed Retraining**: Re-train models using only causally verified predictors (`15b_causal_informed_retraining.R`).
5.  **Conservation Planning**: Generate CATE maps and future suitability projections (`11d_cate_maps.R`, `15_future_env_projection.R`).

## 📊 Key Findings

*   **Efficiency**: Causal selection reduced predictors from **47 to ~29**, improving model transferability.
*   **Stability**: Identified stable causal modules: *Topography → Climate → Land Cover/Soil*.
*   **Uncertainty**: Found that structural (algorithmic) uncertainty significantly exceeds climate scenario uncertainty in mid-century projections.

## 🛠️ Requirements

*   **R** (version >= 4.0.0)
*   **Key R Packages**: `maxnet`, `randomForest`, `mgcv`, `nnet`, `pcalg`, `bnlearn`, `grf`, `DoubleML`, `terra`, `sf`.
*   **Python** (optional, for specific visualizations): `matplotlib`, `seaborn`.

## 📜 Citation

If you use this code or methodology, please cite the associated manuscript:

> *Causal inference reveals mechanism-driven freshwater fish distribution under climate change: a river network perspective.* (In Preparation/Review)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
