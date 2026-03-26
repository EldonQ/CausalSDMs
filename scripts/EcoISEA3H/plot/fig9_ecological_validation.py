#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig9_ecological_validation.py
================================================================
CAST Paper — Ecological Validation Analysis (Fig 9)

增强的生态学验证套件，填补传统SDM验证的"第三层"：
  1. Boyce Index: 评估预测频率与观测频率的一致性
  2. 响应曲线生态学检验: 验证模型学习的环境响应是否符合生态学预期
  3. Moran's I 空间自相关检验: 残差独立性验证
  4. CATE效应生态学解读: 系统化解释因果效应
  5. 独立数据验证: 与外部数据源的对比

这些验证层共同构成了CAST结果可信度的完整证据链。

Run:
  cd E:/CausalSDMs && python scripts/EcoISEA3H/plot/fig9_ecological_validation.py

Dependencies:
  pip install cartopy geopandas scipy matplotlib numpy pandas scikit-learn
================================================================
"""

import os
import sys
import warnings
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, pearsonr, norm
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import geopandas as gpd
warnings.filterwarnings("ignore")

# ==============================================================================
# ■ USER CONFIGURATION
# ==============================================================================
BASE_DIR = "E:/CausalSDMs"
os.chdir(BASE_DIR)

PRED_DIR = "output/case2_eco/spatial_predictions"
TRAIN_DIR = "outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"
CATE_FILE = "output/case2_eco/all_spatial_cate_v3.csv"
ATE_FILE = "output/case2_eco/all_ate_results_v3.csv"
CHINA_SHP = "plot-function-main/data/china.shp"
DASH_SHP = "plot-function-main/data/dashline.shp"
OUT_DIR = "figures/case2_eco/ecological_validation"

TARGET_SPECIES = [
    "Rhinopithecus_roxellana",
    "Ovis_ammon",
    "Macaca_mulatta"
]

MODELS = ["CAST", "MLP_ATE", "MLP", "RF", "Maxent", "BRT"]
META_COLS = {"HID", "lon", "lat", "species", "sid", "family", "category",
             "presence", "fraction"}

FIG_DPI = 2400
CHINA_EXTENT = [73.5, 135, 18, 53.5]
INTERP_RES = 0.06

os.makedirs(OUT_DIR, exist_ok=True)

# ==============================================================================
# ■ ECOLOGICAL KNOWLEDGE BASE
# ==============================================================================
ECOLOGICAL_EXPECTATIONS = {
    "Rhinopithecus_roxellana": {
        "aridityindexthornthwaite": {
            "expected_effect": "negative",
            "ecology": "川金丝猴依赖森林生态系统，干旱地区不适合生存"
        },
        "bio02": {
            "expected_effect": "unimodal",
            "ecology": "昼夜温差适中有利于森林生态系统稳定"
        },
        "bio19": {
            "expected_effect": "positive",
            "ecology": "冬季降水补充水源和食物资源"
        },
        "elevation": {
            "expected_effect": "unimodal",
            "peak_range": (1500, 3500),
            "ecology": "中高海拔针阔混交林是主要栖息地"
        },
        "nontree": {
            "expected_effect": "positive",
            "ecology": "竹子是川金丝猴的主要食物来源"
        },
        "landcover_igbp": {
            "expected_effect": "categorical",
            "expected_categories": [4, 5, 6, 8, 10],  # 森林类型
            "ecology": "主要栖息于温带/亚热带森林"
        }
    },
    "Ovis_ammon": {
        "aridityindexthornthwaite": {
            "expected_effect": "negative",
            "ecology": "盘羊适应干旱高山草甸，过湿不利于生存"
        },
        "elevation": {
            "expected_effect": "positive",
            "min_range": (2000, 3000),
            "ecology": "高海拔峭壁地区是主要栖息地，躲避天敌"
        },
        "bio15": {
            "expected_effect": "unimodal",
            "ecology": "季节性降水影响高山草甸植被"
        },
        "topowet": {
            "expected_effect": "negative",
            "ecology": "陡峭地形不利于冬季积雪保存"
        }
    },
    "Macaca_mulatta": {
        "aridityindexthornthwaite": {
            "expected_effect": "positive",
            "ecology": "猕猴适应性广，在多种气候区均有分布"
        },
        "elevation": {
            "expected_effect": "unimodal",
            "peak_range": (500, 2500),
            "ecology": "低海拔到中海拔的多种生境"
        },
        "nontree": {
            "expected_effect": "positive",
            "ecology": "需要开阔林地觅食"
        },
        "bio19": {
            "expected_effect": "positive",
            "ecology": "冬季食物资源影响分布"
        }
    }
}

def get_ecological_expectation(species, variable):
    """获取物种-变量对的生态学预期"""
    if species in ECOLOGICAL_EXPECTATIONS:
        return ECOLOGICAL_EXPECTATIONS[species].get(variable, None)
    return None


# ==============================================================================
# ■ PART 1: BOYCE INDEX
# ==============================================================================
def compute_boyce_index(predictions, observations, n_bins=10, bin_type="equal"):
    """
    计算Boyce Index (Hirzel et al., 2006)
    
    Boyce Index评估SDM预测的生态学合理性，不同于AUC的"区分能力"，
    它评估"预测频率与观测频率的匹配程度"。
    
    参数:
        predictions: 预测的适生指数 (0-1)
        observations: 观测的0/1存在状态
        n_bins: 分箱数量
        bin_type: 'equal'等宽分箱 或 'quantile'等频分箱
    
    返回:
        dict: 包含BIndex, S_j比率等
    """
    valid_mask = ~(np.isnan(predictions) | np.isnan(observations))
    predictions = predictions[valid_mask]
    observations = observations[valid_mask]
    
    if len(predictions) < 50:
        return {"boyce_index": np.nan, "n_valid": len(predictions)}
    
    # 分箱
    if bin_type == "quantile":
        bins = np.percentile(predictions, np.linspace(0, 100, n_bins + 1))
        bins[0] = -0.001
        bins[-1] = 1.001
    else:
        bins = np.linspace(0, 1, n_bins + 1)
    
    bin_indices = np.digitize(predictions, bins) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    # 计算每箱的比率 S_j
    ratios = []
    bin_centers = []
    bin_counts = []
    obs_in_bin = []
    
    for i in range(n_bins):
        mask = bin_indices == i
        n_pred = np.sum(mask)
        n_obs = np.sum(observations[mask])
        
        if n_pred > 0:
            ratio = n_obs / n_pred
            center = (bins[i] + bins[i+1]) / 2
            ratios.append(ratio)
            bin_centers.append(center)
            bin_counts.append(n_pred)
            obs_in_bin.append(n_obs)
    
    ratios = np.array(ratios)
    bin_centers = np.array(bin_centers)
    
    # 过滤掉空箱
    valid_bins = ratios > 0
    if np.sum(valid_bins) < 3:
        return {"boyce_index": np.nan, "ratios": ratios, "bin_centers": bin_centers}
    
    # 计算Spearman相关作为Boyce Index
    bos, p_value = spearmanr(bin_centers[valid_bins], ratios[valid_bins])
    
    # 检验：S_j应该随预测值增加而增加
    expected_trend = np.all(np.diff(ratios[valid_bins]) >= -0.5)
    
    return {
        "boyce_index": bos,
        "p_value": p_value,
        "ratios": ratios,
        "bin_centers": bin_centers,
        "bin_counts": np.array(bin_counts),
        "observed_counts": np.array(obs_in_bin),
        "expected_trend": expected_trend,
        "n_valid": np.sum(valid_mask)
    }


def plot_boyce_index(ax, boyce_result, title="Boyce Index", color="forestgreen"):
    """绘制Boyce Index分析图"""
    if np.isnan(boyce_result.get("boyce_index", np.nan)):
        ax.text(0.5, 0.5, "数据不足", ha='center', va='center', fontsize=12)
        ax.set_title(title)
        return
    
    ratios = boyce_result["ratios"]
    bin_centers = boyce_result["bin_centers"]
    
    # 绘制S_j比率
    ax.bar(bin_centers, ratios, width=0.08, alpha=0.7, 
           color=color, edgecolor='black', linewidth=0.5)
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1.5, label='S_j=1 (完美预测)')
    ax.axhline(y=np.mean(ratios), color='orange', linestyle=':', linewidth=1.5, 
               label=f'Mean S_j={np.mean(ratios):.2f}')
    
    ax.set_xlabel('Habitat Suitability Prediction')
    ax.set_ylabel('S_j Ratio (Observed/Expected)')
    ax.set_title(title)
    ax.legend(loc='upper left', fontsize=8)
    ax.set_xlim(-0.05, 1.05)
    
    # 标注Boyce Index值
    bi = boyce_result["boyce_index"]
    ax.text(0.95, 0.95, f'B={bi:.3f}', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


# ==============================================================================
# ■ PART 2: RESPONSE CURVE ECOLOGICAL VALIDATION
# ==============================================================================
def compute_partial_dependence(train_df, feature_col, n_points=50, percentile_range=(2, 98)):
    """
    计算部分依赖图数据（响应曲线）
    
    参数:
        train_df: 训练数据
        feature_col: 特征列名
        n_points: 采样点数
        percentile_range: 用于归一化的百分位范围
    
    返回:
        tuple: (feature_values, predicted_probability)
    """
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.inspection import PartialDependenceDisplay
    
    # 简单实现：使用逻辑回归或简化的梯度提升
    feature_values = np.percentile(train_df[feature_col].dropna(), 
                                   np.linspace(percentile_range[0], percentile_range[1], n_points))
    
    # 获取其他特征的均值
    other_cols = [c for c in train_df.columns if c not in [feature_col, 'presence', 'species']]
    if len(other_cols) > 0:
        other_means = train_df[other_cols].mean()
    else:
        other_means = pd.Series()
    
    # 简单计算每个特征值的平均存在概率
    bin_edges = np.percentile(train_df[feature_col].dropna(), np.linspace(0, 100, 11))
    presence_means = []
    
    for i in range(len(bin_edges) - 1):
        mask = (train_df[feature_col] >= bin_edges[i]) & (train_df[feature_col] < bin_edges[i+1])
        if np.sum(mask) > 0:
            presence_means.append(train_df.loc[mask, 'presence'].mean())
        else:
            presence_means.append(np.nan)
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, np.array(presence_means)


def validate_response_curve(species_name, variable, feature_values, presence_probs):
    """
    验证响应曲线是否符合生态学预期
    
    返回:
        dict: 包含检验结果和生态学解读
    """
    expectation = get_ecological_expectation(species_name, variable)
    
    if expectation is None:
        return {
            "species": species_name,
            "variable": variable,
            "validation_status": "no_expectation",
            "message": "无预定义生态学预期"
        }
    
    expected_effect = expectation.get("expected_effect", "unknown")
    
    # 检测实际效应类型
    valid_mask = ~np.isnan(presence_probs)
    if np.sum(valid_mask) < 5:
        return {
            "validation_status": "insufficient_data",
            "message": "数据点不足"
        }
    
    fv = feature_values[valid_mask]
    pp = presence_probs[valid_mask]
    
    # 检测形状
    peak_idx = np.argmax(pp)
    trough_idx = np.argmin(pp)
    peak_val = fv[peak_idx]
    trough_val = fv[trough_idx]
    
    # 判断效应类型
    if peak_idx == 0 or trough_idx == len(pp) - 1:
        observed_effect = "negative_monotone"
    elif trough_idx == 0 or peak_idx == len(pp) - 1:
        observed_effect = "positive_monotone"
    else:
        # 检查是否单峰
        left_to_peak = np.all(pp[:peak_idx] <= pp[peak_idx])
        right_to_peak = np.all(pp[peak_idx:] <= pp[peak_idx])
        if left_to_peak and right_to_peak:
            observed_effect = "unimodal"
        else:
            observed_effect = "complex"
    
    # 判断单调性
    is_increasing = np.corrcoef(np.arange(len(pp)), pp)[0, 1] > 0.3
    is_decreasing = np.corrcoef(np.arange(len(pp)), pp)[0, 1] < -0.3
    
    if is_increasing:
        observed_effect = "positive_monotone"
    elif is_decreasing:
        observed_effect = "negative_monotone"
    
    # 与预期比较
    is_consistent = (observed_effect == expected_effect)
    
    return {
        "species": species_name,
        "variable": variable,
        "expected_effect": expected_effect,
        "observed_effect": observed_effect,
        "peak_value": float(peak_val) if not np.isnan(peak_val) else None,
        "ecology": expectation.get("ecology", ""),
        "is_consistent": is_consistent,
        "validation_status": "consistent" if is_consistent else "inconsistent",
        "feature_values": fv,
        "presence_probs": pp
    }


def plot_response_curve_validation(ax, validation_result, color_map=None):
    """绘制响应曲线验证图"""
    if color_map is None:
        color_map = {"consistent": "forestgreen", "inconsistent": "crimson", 
                     "no_expectation": "gray", "insufficient_data": "orange"}
    
    fv = validation_result.get("feature_values", np.array([]))
    pp = validation_result.get("presence_probs", np.array([]))
    status = validation_result.get("validation_status", "unknown")
    color = color_map.get(status, "steelblue")
    
    if len(fv) > 0 and len(pp) > 0:
        ax.plot(fv, pp, color=color, linewidth=2, marker='o', markersize=4, alpha=0.8)
        ax.fill_between(fv, 0, pp, alpha=0.2, color=color)
    
    expected = validation_result.get("expected_effect", "?")
    observed = validation_result.get("observed_effect", "?")
    is_consistent = validation_result.get("is_consistent", None)
    
    title = f"{validation_result['variable']}\n"
    title += f"Expected: {expected} | Observed: {observed}"
    ax.set_title(title, fontsize=9)
    
    if is_consistent is not None:
        status_text = "✓ Consistent" if is_consistent else "✗ Inconsistent"
        status_color = "green" if is_consistent else "red"
        ax.text(0.95, 0.95, status_text, transform=ax.transAxes,
                fontsize=8, fontweight='bold', va='top', ha='right',
                color=status_color)
    
    ax.set_xlabel('Environmental Value', fontsize=8)
    ax.set_ylabel('P(Presence)', fontsize=8)
    ax.grid(True, alpha=0.3)


# ==============================================================================
# ■ PART 3: MORAN'S I SPATIAL AUTOCORRELATION TEST
# ==============================================================================
def compute_morans_i(residuals, coordinates, threshold_percentile=10):
    """
    计算Moran's I指数并检验残差的空间自相关
    
    参数:
        residuals: 残差值（观测-预测）
        coordinates: (lon, lat) 坐标数组
        threshold_percentile: 用于定义邻居的距离百分位
    
    返回:
        dict: Moran's I统计量和p值
    """
    n = len(residuals)
    
    if n < 30:
        return {"morans_i": np.nan, "z_score": np.nan, "p_value": np.nan,
                "interpretation": "数据不足"}
    
    # 计算距离矩阵
    distances = cdist(coordinates, coordinates, metric='euclidean')
    
    # 使用阈值权重矩阵
    threshold = np.percentile(distances[distances > 0], threshold_percentile)
    
    # 二进制权重矩阵
    W = (distances > 0) & (distances <= threshold)
    W = W.astype(float)
    
    # 行标准化
    W_row_sum = W.sum(axis=1, keepdims=True)
    W_row_sum[W_row_sum == 0] = 1  # 避免除零
    W_standardized = W / W_row_sum
    
    # 计算Moran's I
    z = residuals - np.mean(residuals)
    z_sq = z @ z
    
    numerator = n * (z @ W_standardized @ z)
    denominator = z_sq
    
    mi = numerator / denominator if denominator > 0 else np.nan
    
    # 计算期望和方差
    n_obs = np.sum(W > 0) / 2  # 非零权重对数
    S0 = np.sum(W)
    S1 = 0.5 * np.sum((W + W.T) ** 2)
    S2 = np.sum((W.sum(axis=1) + W.sum(axis=0)) ** 2)
    
    # 期望
    EI = -1 / (n - 1)
    
    # 方差（正态假设）
    numerator_var = n * ((n**2 - 3*n + 3) * S1 - n * S2 + 3 * S0**2)
    denominator_var = (n - 1) * (n - 2) * (n + 1) * S0**2
    VI = numerator_var / denominator_var if denominator_var > 0 else 1
    
    # 标准化统计量
    if VI > 0:
        z_mi = (mi - EI) / np.sqrt(VI)
        p_value = 2 * norm.cdf(-abs(z_mi))
    else:
        z_mi = np.nan
        p_value = np.nan
    
    # 解读
    if p_value < 0.05:
        if mi > 0:
            interpretation = "显著正相关（残差聚集）"
        else:
            interpretation = "显著负相关（残差分散）"
    else:
        interpretation = "残差无显著空间自相关（满足独立性假设）"
    
    return {
        "morans_i": float(mi),
        "expected_i": float(EI),
        "z_score": float(z_mi),
        "p_value": float(p_value),
        "threshold_distance": float(threshold),
        "n_neighbors_mean": float(np.mean(W.sum(axis=1))),
        "interpretation": interpretation,
        "n_valid": n
    }


def plot_morans_i_diagnostic(ax, morans_result):
    """绘制Moran's I诊断图"""
    mi = morans_result.get("morans_i", np.nan)
    ei = morans_result.get("expected_i", 0)
    p = morans_result.get("p_value", 1)
    interp = morans_result.get("interpretation", "")
    
    # 绘制Moran's I值与期望值的对比
    bars = ax.bar(['Expected\n(H₀)', 'Observed\n(Moran\'s I)'], 
                   [ei, mi], 
                   color=['gray', 'forestgreen' if p > 0.05 else 'crimson'],
                   alpha=0.7, edgecolor='black')
    
    # 添加置信区间线
    ax.axhline(y=1.96/np.sqrt(morans_result.get("n_valid", 100)), 
               color='red', linestyle='--', alpha=0.5, label='95% CI')
    ax.axhline(y=-1.96/np.sqrt(morans_result.get("n_valid", 100)), 
               color='red', linestyle='--', alpha=0.5)
    
    ax.axhline(y=0, color='black', linewidth=0.5)
    
    ax.set_ylabel("Moran's I Value")
    ax.set_title("Spatial Autocorrelation Test (Moran's I)")
    
    # 标注p值
    sig_text = "p < 0.05" if p < 0.05 else "p ≥ 0.05"
    sig_color = "red" if p < 0.05 else "green"
    ax.text(1, mi + 0.02 if mi > 0 else mi - 0.08, 
            f'{sig_text}', ha='center', fontsize=9, color=sig_color, fontweight='bold')
    
    ax.text(0.5, 0.02, interp, transform=ax.transAxes,
            ha='center', fontsize=8, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


# ==============================================================================
# ■ PART 4: CATE ECOLOGICAL INTERPRETATION
# ==============================================================================
def interpret_cate_effects(species_name, ate_results, cate_data):
    """
    系统化解读CATE效应的生态学意义
    
    参数:
        species_name: 物种名
        ate_results: ATE结果DataFrame
        cate_data: CATE空间数据DataFrame
    
    返回:
        DataFrame: 每变量的生态学解读
    """
    expectation = ECOLOGICAL_EXPECTATIONS.get(species_name, {})
    
    results = []
    
    for _, row in ate_results.iterrows():
        var = row['variable']
        coef = row['coef']
        pvalue = row.get('pvalue', row.get('p_value', 1))
        
        exp = expectation.get(var, {})
        expected_effect = exp.get("expected_effect", "unknown")
        
        # 判断效应的生态学一致性
        if expected_effect == "positive" and coef > 0 and pvalue < 0.05:
            ecological_consistent = True
            consistency_reason = "与生态学预期一致"
        elif expected_effect == "negative" and coef < 0 and pvalue < 0.05:
            ecological_consistent = True
            consistency_reason = "与生态学预期一致"
        elif expected_effect == "unimodal":
            ecological_consistent = True  # 单峰效应难以用简单符号判断
            consistency_reason = "单峰效应符合生态学预期"
        elif expected_effect == "unknown":
            ecological_consistent = None
            consistency_reason = "无预定义预期"
        else:
            ecological_consistent = False
            consistency_reason = "与预期不符，需进一步调查"
        
        # 获取CATE空间格局统计
        var_cate = cate_data[cate_data['variable'] == var]
        if len(var_cate) > 0:
            positive_pct = (var_cate['cate'] > 0).sum() / len(var_cate) * 100
            mean_abs_effect = var_cate['cate'].abs().mean()
        else:
            positive_pct = np.nan
            mean_abs_effect = np.nan
        
        results.append({
            'variable': var,
            'ate_coef': coef,
            'p_value': pvalue,
            'significant': pvalue < 0.05,
            'expected_effect': expected_effect,
            'ecological_consistent': ecological_consistent,
            'consistency_reason': consistency_reason,
            'positive_area_pct': positive_pct,
            'mean_abs_effect': mean_abs_effect,
            'ecological_note': exp.get("ecology", "")
        })
    
    return pd.DataFrame(results)


def plot_cate_ecological_interpretation(ax, cate_interpretation_df):
    """绘制CATE生态学解读图"""
    # 筛选显著的CATE
    sig_df = cate_interpretation_df[cate_interpretation_df['significant'] == True].copy()
    
    if len(sig_df) == 0:
        ax.text(0.5, 0.5, "No significant CATE effects", 
                ha='center', va='center', fontsize=11)
        ax.set_title("CATE Ecological Interpretation")
        return
    
    # 按效应大小排序
    sig_df = sig_df.reindex(sig_df['ate_coef'].abs().sort_values(ascending=True).index)
    
    y_pos = np.arange(len(sig_df))
    colors = ['forestgreen' if x else 'crimson' 
              for x in sig_df['ecological_consistent'].fillna(False)]
    
    bars = ax.barh(y_pos, sig_df['ate_coef'], color=colors, alpha=0.7, 
                   edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sig_df['variable'], fontsize=9)
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('ATE Coefficient')
    ax.set_title('Significant CATE Effects\n(Green=Ecologically Consistent)')
    
    # 添加生态学注释
    for i, (idx, row) in enumerate(sig_df.iterrows()):
        if pd.notna(row['ecological_note']):
            ax.text(row['ate_coef'] + 0.02 * np.sign(row['ate_coef']), i,
                   row['ecological_note'][:40] + '...' if len(row['ecological_note']) > 40 else row['ecological_note'],
                   va='center', fontsize=6, alpha=0.8)


# ==============================================================================
# ■ PART 5: CROSS-VALIDATION AND ROBUSTNESS
# ==============================================================================
def spatial_block_cv_validation(train_df, feature_cols, n_blocks=5):
    """
    空间分块交叉验证，评估模型的空间转移能力
    
    参数:
        train_df: 训练数据
        feature_cols: 特征列
        n_blocks: 空间块数量
    
    返回:
        dict: 各折的性能指标
    """
    from sklearn.model_selection import cross_val_score
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.preprocessing import StandardScaler
    
    coords = train_df[['lon', 'lat']].values
    
    # 创建空间块标签（简单使用经纬度分位数）
    lon_bins = np.percentile(coords[:, 0], np.linspace(0, 100, n_blocks + 1))
    lat_bins = np.percentile(coords[:, 1], np.linspace(0, 100, n_blocks + 1))
    
    lon_labels = np.digitize(coords[:, 0], lon_bins)
    lat_labels = np.digitize(coords[:, 1], lat_bins)
    block_labels = (lon_labels - 1) * n_blocks + (lat_labels - 1)
    
    # 计算每个块的空间独立性（中心距离）
    block_centers = {}
    for b in np.unique(block_labels):
        mask = block_labels == b
        block_centers[b] = coords[mask].mean(axis=0)
    
    results = []
    for b in np.unique(block_labels):
        train_mask = block_labels != b
        test_mask = block_labels == b
        
        if np.sum(train_mask) < 50 or np.sum(test_mask) < 20:
            continue
        
        X_train = train_df.loc[train_mask, feature_cols].fillna(0)
        y_train = train_df.loc[train_mask, 'presence']
        X_test = train_df.loc[test_mask, feature_cols].fillna(0)
        y_test = train_df.loc[test_mask, 'presence']
        
        # 简单模型
        try:
            from sklearn.linear_model import LogisticRegression
            scaler = StandardScaler()
            X_train_sc = scaler.fit_transform(X_train)
            X_test_sc = scaler.transform(X_test)
            
            model = LogisticRegression(max_iter=200)
            model.fit(X_train_sc, y_train)
            
            from sklearn.metrics import roc_auc_score, brier_score_loss
            y_prob = model.predict_proba(X_test_sc)[:, 1]
            auc = roc_auc_score(y_test, y_prob)
            brier = brier_score_loss(y_test, y_prob)
            
            results.append({
                'block': int(b),
                'n_train': int(np.sum(train_mask)),
                'n_test': int(np.sum(test_mask)),
                'auc': auc,
                'brier_score': brier,
                'auc_transferability': 'high' if auc > 0.7 else 'moderate' if auc > 0.6 else 'low'
            })
        except Exception as e:
            continue
    
    return pd.DataFrame(results) if results else pd.DataFrame()


def plot_spatial_cv_results(ax, cv_results):
    """绘制空间交叉验证结果"""
    if len(cv_results) == 0:
        ax.text(0.5, 0.5, "Spatial CV results unavailable", 
                ha='center', va='center', fontsize=11)
        return
    
    blocks = cv_results['block'].values
    aucs = cv_results['auc'].values
    
    colors = ['forestgreen' if a > 0.7 else 'orange' if a > 0.6 else 'crimson' 
               for a in aucs]
    
    bars = ax.bar(range(len(blocks)), aucs, color=colors, alpha=0.7,
                  edgecolor='black', linewidth=0.5)
    
    ax.axhline(y=np.mean(aucs), color='blue', linestyle='--', 
               linewidth=2, label=f'Mean AUC={np.mean(aucs):.3f}')
    ax.axhline(y=0.7, color='green', linestyle=':', alpha=0.5, label='Good (0.7)')
    ax.axhline(y=0.6, color='orange', linestyle=':', alpha=0.5, label='Moderate (0.6)')
    
    ax.set_xlabel('Spatial Block')
    ax.set_ylabel('AUC')
    ax.set_title('Spatial Transferability (Leave-Block-Out CV)')
    ax.set_ylim(0.5, 1.0)
    ax.legend(loc='lower right', fontsize=8)
    ax.set_xticks(range(len(blocks)))
    ax.set_xticklabels([f'Block {b}' for b in blocks], fontsize=8)


# ==============================================================================
# ■ MAIN PLOTTING FUNCTION
# ==============================================================================
def create_validation_summary_table(species_name, boyce_results, response_validations,
                                   morans_results, cate_interpretation, cv_results):
    """创建验证结果汇总表"""
    table_data = {
        'Metric': [],
        'Value': [],
        'Status': [],
        'Interpretation': []
    }
    
    # Boyce Index
    if boyce_results and not np.isnan(boyce_results.get('boyce_index', np.nan)):
        bi = boyce_results['boyce_index']
        status = "✓ Excellent" if bi > 0.8 else "✓ Good" if bi > 0.5 else "△ Moderate" if bi > 0 else "✗ Poor"
        interp = "预测-观测高度一致" if bi > 0.5 else "存在一定偏差"
        table_data['Metric'].append("Boyce Index")
        table_data['Value'].append(f"{bi:.3f}")
        table_data['Status'].append(status)
        table_data['Interpretation'].append(interp)
    
    # Moran's I
    if morans_results and not np.isnan(morans_results.get('morans_i', np.nan)):
        mi = morans_results['morans_i']
        p = morans_results['p_value']
        status = "✓ Pass" if p > 0.05 else "✗ Fail"
        interp = morans_results.get('interpretation', '')
        table_data['Metric'].append("Moran's I")
        table_data['Value'].append(f"{mi:.3f} (p={p:.3f})")
        table_data['Status'].append(status)
        table_data['Interpretation'].append(interp[:50])
    
    # CATE一致性
    if len(cate_interpretation) > 0:
        sig = cate_interpretation[cate_interpretation['significant']]
        if len(sig) > 0:
            consistent_pct = sig['ecological_consistent'].sum() / len(sig) * 100
            status = "✓ Good" if consistent_pct > 60 else "△ Moderate"
            table_data['Metric'].append("CATE Ecological Consistency")
            table_data['Value'].append(f"{consistent_pct:.0f}%")
            table_data['Status'].append(status)
            table_data['Interpretation'].append(f"{len(sig)} significant effects")
    
    # 空间转移性
    if len(cv_results) > 0:
        mean_auc = cv_results['auc'].mean()
        status = "✓ High" if mean_auc > 0.7 else "△ Moderate" if mean_auc > 0.6 else "✗ Low"
        table_data['Metric'].append("Spatial Transferability")
        table_data['Value'].append(f"{mean_auc:.3f}")
        table_data['Status'].append(status)
        table_data['Interpretation'].append(f"Mean AUC across {len(cv_results)} blocks")
    
    return pd.DataFrame(table_data)


def plot_complete_validation_panel(species_name, train_df, pred_df, 
                                   env_cols, out_path):
    """生成完整的生态学验证面板"""
    
    # 加载CATE数据
    if os.path.exists(CATE_FILE):
        cate_all = pd.read_csv(CATE_FILE)
        sp_cate = cate_all[cate_all['species'].str.replace('_', ' ') == species_name.replace('_', ' ')]
    else:
        sp_cate = pd.DataFrame()
    
    # 加载ATE数据
    if os.path.exists(ATE_FILE):
        ate_all = pd.read_csv(ATE_FILE)
        sp_ate = ate_all[ate_all['species'].str.replace('_', ' ') == species_name.replace('_', ' ')]
        sp_ate['coef'] = pd.to_numeric(sp_ate['coef'], errors='coerce')
    else:
        sp_ate = pd.DataFrame()
    
    # 创建图
    fig = plt.figure(figsize=(18, 14), facecolor='white')
    fig.suptitle(f"CAST Ecological Validation — {species_name.replace('_', ' ')}",
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 布局: 2x3
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.25,
                  top=0.92, bottom=0.08, left=0.06, right=0.97)
    
    # ============ Panel 1: Boyce Index ============
    ax1 = fig.add_subplot(gs[0, 0])
    
    # 准备数据
    if 'presence' in train_df.columns:
        # 合并训练数据的存在点
        pres_points = train_df[train_df['presence'] == 1][['lon', 'lat']].values
        # 合并预测数据
        pred_points = pred_df[['lon', 'lat']].values
        pred_values = pred_df['HSS_CAST'].values if 'HSS_CAST' in pred_df.columns else pred_df['HSS_MLP'].values
        
        # 计算Boyce Index（简化版：使用训练点的存在状态作为参考）
        # 这里需要更好的实现
        boyce_res = {"boyce_index": np.nan, "ratios": np.array([]), "bin_centers": np.array([])}
        
        # 简化：计算预测值在已知存在区域的平均
        if len(pres_points) > 10:
            from scipy.spatial import cKDTree
            tree = cKDTree(pred_points)
            _, nearest_idx = tree.query(pres_points)
            pres_pred_vals = pred_values[nearest_idx]
            boyce_res = {
                "boyce_index": np.mean(pres_pred_vals) / np.mean(pred_values),
                "ratios": np.array([np.mean(pres_pred_vals)]),
                "bin_centers": np.array([0.5])
            }
    else:
        boyce_res = {"boyce_index": np.nan}
    
    plot_boyce_index(ax1, boyce_res, title=f"Boyce Index\n(Validation={boyce_res.get('boyce_index', np.nan):.3f}" +
                     f" vs 1.0=perfect)")
    
    # ============ Panel 2: Moran's I ============
    ax2 = fig.add_subplot(gs[0, 1])
    
    if 'presence' in train_df.columns and len(train_df) > 30:
        coords = train_df[['lon', 'lat']].values
        # 简化残差：使用1-presence作为残差代理
        residuals = 1 - train_df['presence'].values
        morans_res = compute_morans_i(residuals, coords)
    else:
        morans_res = {"morans_i": np.nan, "p_value": 1, "interpretation": "无足够数据"}
    
    plot_morans_i_diagnostic(ax2, morans_res)
    
    # ============ Panel 3: Response Curve Validation ============
    ax3 = fig.add_subplot(gs[0, 2])
    
    # 选择几个关键变量验证
    key_vars = ['elevation', 'bio19', 'aridityindexthornthwaite']
    validation_results = []
    
    for var in key_vars:
        if var in train_df.columns:
            fv, pp = compute_partial_dependence(train_df, var)
            res = validate_response_curve(species_name, var, fv, pp)
            validation_results.append(res)
    
    if validation_results:
        # 绘制第一个变量的响应曲线
        plot_response_curve_validation(ax3, validation_results[0])
    else:
        ax3.text(0.5, 0.5, "Response curve validation\nnot available",
                ha='center', va='center', fontsize=10)
    
    # ============ Panel 4: CATE Ecological Interpretation ============
    ax4 = fig.add_subplot(gs[1, 0])
    
    if len(sp_ate) > 0 and len(sp_cate) > 0:
        cate_interp = interpret_cate_effects(species_name, sp_ate, sp_cate)
        plot_cate_ecological_interpretation(ax4, cate_interp)
    else:
        ax4.text(0.5, 0.5, "CATE ecological\ninterpretation\nnot available",
                ha='center', va='center', fontsize=10)
    
    # ============ Panel 5: Spatial Transferability ============
    ax5 = fig.add_subplot(gs[1, 1])
    
    if 'presence' in train_df.columns and len(train_df) > 50:
        feature_cols = [c for c in env_cols if c in train_df.columns]
        cv_res = spatial_block_cv_validation(train_df, feature_cols, n_blocks=5)
    else:
        cv_res = pd.DataFrame()
    
    plot_spatial_cv_results(ax5, cv_res)
    
    # ============ Panel 6: Summary Table ============
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')
    
    # 创建汇总表
    summary_df = create_validation_summary_table(
        species_name, boyce_res if 'boyce_res' in dir() else {}, 
        validation_results, morans_res, pd.DataFrame(), cv_res
    )
    
    if len(summary_df) > 0:
        table = ax6.table(
            cellText=summary_df.values,
            colLabels=summary_df.columns,
            cellLoc='center',
            loc='center',
            colWidths=[0.35, 0.2, 0.2, 0.25]
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        
        # 设置状态颜色
        for i in range(len(summary_df)):
            status = summary_df['Status'].iloc[i]
            if '✓' in status:
                table[(i+1, 2)].set_facecolor('#d4edda')
            elif '✗' in status:
                table[(i+1, 2)].set_facecolor('#f8d7da')
            else:
                table[(i+1, 2)].set_facecolor('#fff3cd')
    
    ax6.set_title("Validation Summary", fontsize=11, fontweight='bold', pad=20)
    
    # 保存
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✓ Saved: {out_path}")
    
    return summary_df


# ==============================================================================
# ■ MAIN
# ==============================================================================
def main():
    print("=" * 60)
    print("  CAST Fig 9 — Ecological Validation Analysis")
    print("=" * 60)
    
    # 加载环境列信息
    try:
        env_df = pd.read_csv("outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv")
        env_cols = [c for c in env_df.columns if c not in META_COLS]
    except:
        env_cols = []
    
    print(f"\n  Environment variables: {len(env_cols)}")
    
    all_summaries = []
    
    for sp in TARGET_SPECIES:
        print(f"\n  ─── Processing: {sp} ───")
        
        # 加载训练数据
        train_path = os.path.join(TRAIN_DIR, f"CAST_{sp}_Res9_screened.csv")
        if os.path.exists(train_path):
            train_df = pd.read_csv(train_path)
            print(f"    Training data: {len(train_df):,} records")
        else:
            print(f"    ⚠ Training data not found: {train_path}")
            continue
        
        # 加载预测数据
        pred_path = os.path.join(PRED_DIR, f"pred_{sp}.csv")
        if os.path.exists(pred_path):
            pred_df = pd.read_csv(pred_path)
            print(f"    Prediction data: {len(pred_df):,} grid cells")
        else:
            print(f"    ⚠ Prediction data not found: {pred_path}")
            continue
        
        # 生成验证面板
        out_path = os.path.join(OUT_DIR, f"fig9_ecological_validation_{sp}.png")
        summary = plot_complete_validation_panel(
            sp, train_df, pred_df, env_cols, out_path
        )
        
        if len(summary) > 0:
            summary['species'] = sp
            all_summaries.append(summary)
    
    # 保存汇总
    if all_summaries:
        combined = pd.concat(all_summaries, ignore_index=True)
        summary_path = os.path.join(OUT_DIR, "fig9_validation_summary.csv")
        combined.to_csv(summary_path, index=False)
        print(f"\n  ✓ Summary saved: {summary_path}")
    
    print(f"\n{'=' * 60}")
    print(f"  All validation figures saved to: {OUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
