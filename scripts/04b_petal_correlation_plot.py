#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名称: 04b_petal_correlation_plot.py
功能说明: 花瓣状相关性热图 - 展示环境变量分组间的相关性
          (基于 variables_selected_47.csv 的元数据进行分组和标记)
输入文件: 
  - output/04_collinearity/collinearity_removed.csv (数据)
  - scripts/variables_selected_47.csv (元数据：分组、标签)
输出文件: figures/04_collinearity/petal_correlation_plot.png/svg
作者: Nature级别科研项目
日期: 2025-12-03
==============================================================================
"""

# =============================================================================
# 1. 库的导入
# =============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Wedge, Circle
from matplotlib.lines import Line2D
from scipy import stats
import matplotlib
import os

# 设置字体与SVG格式
matplotlib.rcParams['svg.fonttype'] = 'none' # 文本作为文本保存，便于编辑
plt.rcParams['font.family'] = 'Arial'  # Nature期刊要求Arial

# =============================================================================
# 2. 颜色库设置
# =============================================================================
# 针对4个特定分组的配色方案
COLOR_THEMES = {
    'default': {
        'group_colors': {
            'Topography & Hydrology': '#984EA3',  # 紫色系
            'Upstream Bioclim': '#E41A1C',        # 红色系
            'Land Cover': '#4DAF4A',              # 绿色系
            'Soil': '#A65628'                     # 棕色系
        },
        'heatmap_colors': ['#2166AC', '#F7F7F7', '#B2182B'],  # 蓝-白-红 (冷暖色调)
        'group_label_color': 'black'
    }
}

# =============================================================================
# 3. 绘图前的准备
# =============================================================================
# 选择分析方法 (spearman, pearson, kendall)
selected_method = 'pearson' # 连续变量通常用Pearson，也可选Spearman

# 路径设置
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path = os.path.join(base_dir, "output", "04_collinearity", "collinearity_removed.csv")
metadata_path = os.path.join(base_dir, "scripts", "variables_selected_47.csv")
output_directory = os.path.join(base_dir, "figures", "04_collinearity")

# 确保输出目录存在
os.makedirs(output_directory, exist_ok=True)

# 提取配色方案
select_color = COLOR_THEMES['default']

print("=" * 80)
print("花瓣状相关性热图绘制 - 4组变量相关性分析")
print("=" * 80)
print(f"分析方法: {selected_method.upper()}")
print("")

# =============================================================================
# 4. 数据读取与变量分组
# =============================================================================
print("步骤 1/6: 读取数据与元数据...")

# 1. 读取元数据 (获取分组和标签)
if not os.path.exists(metadata_path):
    raise FileNotFoundError(f"未找到元数据文件: {metadata_path}")

meta_df = pd.read_csv(metadata_path)
# 创建映射字典
var_to_group = dict(zip(meta_df['variable'], meta_df['heatmap_group']))
var_to_label = dict(zip(meta_df['variable'], meta_df['plot_label']))

# 获取有序的分组列表 (确保顺序固定)
# 期望顺序: Topography & Hydrology, Upstream Bioclim, Land Cover, Soil
desired_order = ['Topography & Hydrology', 'Upstream Bioclim', 'Land Cover', 'Soil']
unique_groups = [g for g in desired_order if g in meta_df['heatmap_group'].unique()]

print(f"  - 元数据加载成功，包含 {len(meta_df)} 个变量")
print(f"  - 变量分组: {', '.join(unique_groups)}")

# 2. 读取实际数据
if not os.path.exists(data_path):
    raise FileNotFoundError(f"未找到数据文件: {data_path}")

data = pd.read_csv(data_path)
print(f"  - 数据加载成功，样本数: {len(data)}")

# 3. 构建变量分组字典 (仅包含数据中存在的变量)
data_cols = data.columns.tolist()
var_groups = {}

for group in unique_groups:
    # 找出该组的所有变量
    group_vars_meta = meta_df[meta_df['heatmap_group'] == group]['variable'].tolist()
    # 筛选出数据中实际存在的变量
    existing_vars = [v for v in group_vars_meta if v in data_cols]
    
    if existing_vars:
        var_groups[group] = existing_vars

# 统计每组变量数
print("\n变量分组统计:")
for group_name, group_vars in var_groups.items():
    print(f"  - {group_name}: {len(group_vars)} 个变量")

# =============================================================================
# 5. 计算分组内变量的汇总代表值（用于跨组相关分析）
# =============================================================================
print("\n步骤 2/6: 计算各组代表值（标准化后均值）...")

group_representatives = {}
env_data = data[sum(var_groups.values(), [])] # 仅提取相关环境变量

for group_name, group_vars in var_groups.items():
    if len(group_vars) > 0:
        # 标准化后取均值作为该组的代表值 (PC1也是一种选择，但均值更直观)
        subset = env_data[group_vars]
        group_data_scaled = (subset - subset.mean()) / subset.std()
        group_representatives[group_name] = group_data_scaled.mean(axis=1)

group_data_df = pd.DataFrame(group_representatives)

# =============================================================================
# 6. 计算相关性矩阵
# =============================================================================
print(f"\n步骤 3/6: 计算{selected_method.upper()}相关系数...")

# 每个分组内变量 vs. 所有分组的相关性
all_correlation_data = {}
group_names = list(var_groups.keys())

for group_name in group_names:
    features = var_groups[group_name]
    if len(features) == 0:
        continue
    
    # 目标：其他分组的代表值
    targets = [g for g in group_names if g != group_name]
    
    # 计算相关性矩阵
    n_features = len(features)
    n_targets = len(targets)
    
    correlation_matrix = np.zeros((n_features, n_targets))
    p_value_matrix = np.zeros((n_features, n_targets))
    
    for i, feature_name in enumerate(features):
        for j, target_name in enumerate(targets):
            feature_col = pd.to_numeric(env_data[feature_name], errors='coerce')
            target_col = pd.to_numeric(group_data_df[target_name], errors='coerce')
            
            combined = pd.concat([feature_col, target_col], axis=1).dropna()
            
            if len(combined) < 2:
                corr, p_value = np.nan, np.nan
            else:
                if selected_method == 'spearman':
                    corr, p_value = stats.spearmanr(combined.iloc[:, 0], combined.iloc[:, 1])
                elif selected_method == 'pearson':
                    corr, p_value = stats.pearsonr(combined.iloc[:, 0], combined.iloc[:, 1])
                else:
                    corr, p_value = stats.kendalltau(combined.iloc[:, 0], combined.iloc[:, 1])
            
            correlation_matrix[i, j] = corr
            p_value_matrix[i, j] = p_value
    
    # 保存为DataFrame
    df_corr = pd.DataFrame(correlation_matrix, index=features, columns=targets)
    df_sig = pd.DataFrame(p_value_matrix < 0.05, index=features, columns=targets)
    
    all_correlation_data[group_name] = {
        'correlation_df': df_corr,
        'p_value_df': df_sig
    }
    
    print(f"  - {group_name}: {n_features} 个变量 vs. {n_targets} 个目标组")

# =============================================================================
# 7. 绘图函数
# =============================================================================
def create_full_ring_plot(all_data, all_feature_names, all_target_names, 
                         color_palette, sector_params, var_labels):
    """
    创建花瓣状相关性热图
    """
    
    # 创建画布
    fig, ax = plt.subplots(figsize=(20, 20), subplot_kw={'aspect': 'equal'})
    ax.axis('off')
    
    # 设置颜色映射
    heatmap_colors_value = color_palette['heatmap_colors']
    cmap = LinearSegmentedColormap.from_list(
        "custom_cmap", 
        list(zip([0.0, 0.5, 1.0], heatmap_colors_value))
    )
    norm = plt.Normalize(vmin=-1, vmax=1)
    
    # 获取分组名称
    group_names_list = list(all_data.keys())
    group_legend_colors = [color_palette['group_colors'].get(g, '#333333') for g in group_names_list]
    
    # 遍历每个分组绘制扇区
    for idx, group_name in enumerate(group_names_list):
        features = all_feature_names[group_name]
        df = all_data[group_name]['correlation_df']
        df_sig = all_data[group_name]['p_value_df']
        current_targets = all_target_names[group_name]
        
        # 扇区角度
        start_angle_deg = sector_params[group_name]['start']
        end_angle_deg = sector_params[group_name]['end']
        
        # 特征变量的角度分布
        theta_deg = np.linspace(start_angle_deg, end_angle_deg, len(features))
        theta_rad = np.deg2rad(theta_deg)
        angle_span_deg = abs(end_angle_deg - start_angle_deg) / len(features) * 0.92
        
        current_group_color = group_legend_colors[idx]
        
        # 定义半径 (从内向外)
        radii = np.arange(3, 3 + len(current_targets))
        
        # 绘制每一层 (每个目标组一圈)
        for i, target_name in enumerate(current_targets):
            r_inner = radii[i]
            r_outer = radii[i] + 0.95
            
            values = df[target_name]
            sig_values = df_sig[target_name]
            cell_colors = cmap(norm(values))
            
            # 绘制每个小方格（扇形）
            for j in range(len(features)):
                theta_start = theta_deg[j] - angle_span_deg / 2
                theta_end = theta_deg[j] + angle_span_deg / 2
                
                wedge = Wedge(
                    center=(0, 0),
                    r=r_outer,
                    theta1=theta_start,
                    theta2=theta_end,
                    width=0.92,
                    facecolor=cell_colors[j],
                    edgecolor='#E6E6E6',
                    linewidth=0.6
                )
                ax.add_patch(wedge)
                
                # 添加相关系数文本 (仅在显著或相关性强时显示，避免杂乱)
                # 这里为了清晰，仅显示显著的 * 号，或者只显示强相关的数值
                # 用户要求高DPI，可以显示数值
                text_angle_rad = theta_rad[j]
                text_radius = r_inner + 0.45
                x = text_radius * np.cos(text_angle_rad)
                y = text_radius * np.sin(text_angle_rad)
                
                val = values.iloc[j]
                if not np.isnan(val):
                    sig_marker = '*' if sig_values.iloc[j] else ''
                    # 仅显示数值，若显著加粗或加星
                    text_val = f'{val:.2f}'
                    
                    # 调整文字旋转
                    rot = theta_deg[j] - 90 if np.cos(text_angle_rad) > -0.01 else theta_deg[j] + 90
                    
                    # 字体大小根据变量多少调整
                    font_size = 5 if len(features) > 15 else 6
                    
                    ax.text(
                        x, y, text_val + sig_marker,
                        ha='center', va='center',
                        fontsize=font_size, rotation=rot,
                        color='white' if abs(val) > 0.6 else 'black'
                    )
        
        # 组内环形分隔线
        for r in radii:
            circ = Circle((0, 0), r + 0.95, fill=False, edgecolor='#EDEDED', linewidth=0.5, alpha=0.9)
            ax.add_patch(circ)
        
        # 添加特征变量标签（最外圈）
        label_radius = radii.max() + 0.5
        for i in range(len(features)):
            text_angle_rad = theta_rad[i]
            # 标签向外延伸一点
            x = (label_radius + 0.5) * np.cos(text_angle_rad)
            y = (label_radius + 0.5) * np.sin(text_angle_rad)
            rot = theta_deg[i] if np.cos(text_angle_rad) > -0.01 else theta_deg[i] + 180
            
            # 使用 plot_label
            raw_name = features[i]
            label_text = var_labels.get(raw_name, raw_name)
            
            ax.text(
                x, y, label_text,
                ha='left' if np.cos(text_angle_rad) > -0.01 else 'right', 
                va='center',
                fontsize=7, rotation=rot,
                color=current_group_color,
                fontweight='bold'
            )
        
        # 添加分组标签 (大标题)
        group_label_angle_deg = (start_angle_deg + end_angle_deg) / 2
        group_label_angle_rad = np.deg2rad(group_label_angle_deg)
        # 放在更外层
        group_label_radius = radii.max() + 4.5
        x = group_label_radius * np.cos(group_label_angle_rad)
        y = group_label_radius * np.sin(group_label_angle_rad)
        
        ax.text(
            x, y, group_name,
            ha='center', va='center',
            fontsize=16, fontweight='bold',
            color=current_group_color,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', 
                     edgecolor=current_group_color, linewidth=1.5, alpha=0.9)
        )

        # 分组径向分隔线（起止角）
        for ang in (start_angle_deg, end_angle_deg):
            a = np.deg2rad(ang)
            x_end = (group_label_radius + 2) * np.cos(a)
            y_end = (group_label_radius + 2) * np.sin(a)
            ax.plot([0, x_end], [0, y_end], color='#DDDDDD', linewidth=1.0, zorder=0)
    
    # 创建图例（目标分组）
    # 由于每个组的目标不同，这里做一个统一的图例说明形状代表什么不太容易
    # 但我们的逻辑是：每一圈代表一个目标组。
    # 可以在图的角落列出每个组对应的颜色
    
    legend_elements = []
    for g_name in group_names_list:
        color = color_palette['group_colors'].get(g_name, 'black')
        legend_elements.append(Line2D([0], [0], marker='o', color='w', label=g_name,
                              markerfacecolor=color, markersize=10))
    
    ax.legend(handles=legend_elements, loc='center', title="Variable Groups", 
              fontsize=10, title_fontsize=12, frameon=False)
    
    # 添加颜色条
    cax = fig.add_axes([0.3, 0.05, 0.4, 0.015])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(f"Correlation Coefficient", size=14)
    cbar.ax.tick_params(labelsize=10)
    
    # 设置坐标轴范围
    max_r = max([len(all_target_names[g]) for g in group_names_list]) + 3 + 6
    ax.set_xlim(-max_r, max_r)
    ax.set_ylim(-max_r, max_r)
    
    plt.tight_layout()
    return fig

# =============================================================================
# 8. 执行绘图
# =============================================================================
print("\n步骤 4/6: 绘制花瓣状热图...")

# 设置扇区参数（360度均分）
n_groups = len(group_names)
angle_per_group = 360 / n_groups
gap = 5  # 分组间隙

sector_params = {}
for i, group_name in enumerate(group_names):
    # 调整起始角度，使第一个组在右上方
    start_angle = i * angle_per_group + gap/2 
    end_angle = (i + 1) * angle_per_group - gap/2
    
    sector_params[group_name] = {
        'start': start_angle,
        'end': end_angle
    }

# 准备目标字典（每组对应其他组）
all_target_names = {}
for group_name in group_names:
    # 目标组的顺序可以固定，也可以动态
    all_target_names[group_name] = [g for g in group_names if g != group_name]

# 绘制
fig = create_full_ring_plot(
    all_data=all_correlation_data,
    all_feature_names=var_groups,
    all_target_names=all_target_names,
    color_palette=select_color,
    sector_params=sector_params,
    var_labels=var_to_label
)

# =============================================================================
# 9. 保存结果
# =============================================================================
print("\n步骤 5/6: 保存图表...")

# 保存PNG（高分辨率 2400 DPI）
png_path = os.path.join(output_directory, "petal_correlation_plot.png")
fig.savefig(png_path, dpi=2400, bbox_inches='tight', facecolor='white')
print(f"  - 已保存: {png_path}")

# 保存SVG（矢量图）
svg_path = os.path.join(output_directory, "petal_correlation_plot.svg")
fig.savefig(svg_path, bbox_inches='tight', facecolor='white')
print(f"  - 已保存: {svg_path}")

# 额外导出：将每个变量组的相关性矩阵与显著性矩阵保存为CSV
csv_dir = os.path.join(output_directory, "petal_tables")
os.makedirs(csv_dir, exist_ok=True)
for group_name in group_names:
    corr_df = all_correlation_data[group_name]['correlation_df']
    sig_df = all_correlation_data[group_name]['p_value_df']
    
    # 清理文件名中的特殊字符
    safe_name = group_name.replace(" ", "_").replace("&", "and")
    
    corr_csv = os.path.join(csv_dir, f"{safe_name}_correlation.csv")
    sig_csv = os.path.join(csv_dir, f"{safe_name}_significance.csv")
    corr_df.to_csv(corr_csv, index=True)
    sig_df.to_csv(sig_csv, index=True)
print(f"  - 已保存变量组相关性与显著性表: {csv_dir}")

plt.close()

print("\n步骤 6/6: 总结输出...")
print("\n" + "=" * 80)
print("花瓣状相关性热图绘制完成!")
print("=" * 80)
print(f"\n变量组数量: {len(group_names)}")
print(f"变量组: {', '.join(group_names)}")
print(f"\n输出文件:")
print(f"  - PNG: {png_path}")
print(f"  - SVG: {svg_path}")
print(f"  - CSV表: {csv_dir}")

print("\n版式: Arial 字体, 4组配色, 2400 dpi")
