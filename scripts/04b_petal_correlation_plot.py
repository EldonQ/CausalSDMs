#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名称: 04b_petal_correlation_plot.py
功能说明: 花瓣状相关性热图 - 展示环境变量分组间的相关性
参考: 公众号花瓣状热图教程
输入文件: ../output/04_collinearity/collinearity_removed.csv
输出文件: ../figures/04_collinearity/petal_correlation.png/pdf
作者: Nature级别科研项目
日期: 2025-10-14
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

# 设置字体与PDF格式
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'  # Nature期刊要求Arial

# =============================================================================
# 2. 颜色库设置
# =============================================================================
COLOR_THEMES = {
    1: {
        'group_colors': {
            'Temperature': '#E41A1C',
            'Precipitation': '#377EB8',
            'Hydroclimatic': '#4DAF4A',
            'Topography': '#984EA3',
            'LandCover': '#FF7F00',
            'Soil': '#A65628',
            'Geology': '#F781BF'
        },
        'heatmap_colors': ['#2166AC', '#F7F7F7', '#B2182B'],  # 蓝-白-红
        'group_label_color': 'white'
    },
    2: {
        'group_colors': {
            'Temperature': '#D73027',
            'Precipitation': '#4575B4',
            'Hydroclimatic': '#91BFDB',
            'Topography': '#FC8D59',
            'LandCover': '#FEE090',
            'Soil': '#E0F3F8',
            'Geology': '#FFFFBF'
        },
        'heatmap_colors': 'RdYlBu_r',
        'group_label_color': 'black'
    },
    3: {
        'group_colors': {
            'Temperature': '#D73027',
            'Precipitation': '#1f77b4',
            'Hydroclimatic': '#2ca02c',
            'Topography': '#9467bd',
            'LandCover': '#ff7f0e',
            'Soil': '#8c564b',
            'Geology': '#e377c2'
        },
        'heatmap_colors': ['#2b6cb0', '#f2f2f2', '#f2c94c'],
        'group_label_color': 'black'
    }
}

# =============================================================================
# 3. 绘图前的准备
# =============================================================================
# 选择配色方案（3 更接近参考图）
selected_scheme = 3

# 选择分析方法 (spearman, pearson, kendall)
selected_method = 'spearman'

# 输入文件的地址，输出结果的路径
data_directory = r"E:\SDM01\output\04_collinearity"
output_directory = r"E:\SDM01\figures\04_collinearity"

# 确保输出目录存在
os.makedirs(output_directory, exist_ok=True)

# 从颜色库里提取配色方案
select_color = COLOR_THEMES.get(selected_scheme, COLOR_THEMES[1])

print("=" * 80)
print("花瓣状相关性热图绘制 - 变量组相关性分析")
print("=" * 80)
print(f"分析方法: {selected_method.upper()}")
print(f"配色方案: {selected_scheme}")
print("")

# =============================================================================
# 4. 数据读取与变量分组
# =============================================================================
print("步骤 1/6: 读取数据并按类型进行变量分组...")

# 读取数据
data_path = os.path.join(data_directory, "collinearity_removed.csv")
data = pd.read_csv(data_path)

# 提取环境变量（排除前5列的id, species, lon, lat, source和最后的presence列）
all_cols = data.columns.tolist()
# 找到presence列的位置
presence_idx = [i for i, c in enumerate(all_cols) if 'presence' in c.lower()]
if len(presence_idx) > 0:
    env_columns = all_cols[5:presence_idx[0]]
else:
    env_columns = all_cols[5:]

env_data = data[env_columns].apply(pd.to_numeric, errors='coerce')

print(f"  - 总变量数: {len(env_columns)}")
print(f"  - 样本数: {len(data)}")

# 定义变量分组（按真实筛选后的变量）
var_groups = {
    'Temperature': [col for col in env_columns if col.startswith('tmin_avg') or col.startswith('tmax_avg')],
    'Precipitation': [col for col in env_columns if col.startswith('prec_sum')],
    'Hydroclimatic': [col for col in env_columns if col.startswith('hydro_avg')],
    'Topography': [col for col in env_columns if col.startswith(('dem_avg', 'slope_avg', 'flow_'))],
    'LandCover': [col for col in env_columns if col.startswith('lc_avg')],
    'Soil': [col for col in env_columns if col.startswith('soil_avg')],
    'Geology': [col for col in env_columns if col.startswith('geo_wsum')]
}

# 移除空组
var_groups = {k: v for k, v in var_groups.items() if len(v) > 0}

# 统计每组变量数
print("\n变量分组统计:")
for group_name, group_vars in var_groups.items():
    print(f"  - {group_name}: {len(group_vars)} 个变量")

# =============================================================================
# 5. 计算分组内变量的汇总代表值（用于跨组相关分析）
# =============================================================================
print("\n步骤 2/6: 计算各组代表值（标准化后均值）...")

group_representatives = {}
for group_name, group_vars in var_groups.items():
    if len(group_vars) > 0:
        # 标准化后取均值作为该组的代表值
        group_data_scaled = (env_data[group_vars] - env_data[group_vars].mean()) / env_data[group_vars].std()
        group_representatives[group_name] = group_data_scaled.mean(axis=1)

group_data_df = pd.DataFrame(group_representatives)

# =============================================================================
# 6. 计算相关性矩阵
# =============================================================================
print(f"\n步骤 3/6: 计算{selected_method.upper()}相关系数...")

# 每个分组内变量 vs. 所有分组的相关性
# 这里我们计算：每组的各个变量 vs. 其他组的代表值
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
                         color_palette, sector_params):
    """
    创建花瓣状相关性热图
    
    参数:
    - all_data: 所有分组的相关性数据
    - all_feature_names: 各分组的特征变量名
    - all_target_names: 各分组的目标变量名
    - color_palette: 颜色方案
    - sector_params: 扇区参数（起始角度等）
    """
    
    # 创建画布
    fig, ax = plt.subplots(figsize=(24, 24), subplot_kw={'aspect': 'equal'})
    ax.axis('off')
    
    # 设置颜色映射
    heatmap_colors_value = color_palette['heatmap_colors']
    if isinstance(heatmap_colors_value, str):
        cmap = plt.get_cmap(heatmap_colors_value)
    else:
        cmap = LinearSegmentedColormap.from_list(
            "custom_cmap", 
            list(zip([0.0, 0.5, 1.0], heatmap_colors_value))
        )
    norm = plt.Normalize(vmin=-1, vmax=1)
    
    # 获取分组名称
    group_names_list = list(all_data.keys())
    group_legend_colors = [color_palette['group_colors'][g] for g in group_names_list]
    
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
        
        # 定义半径
        radii = np.arange(3, 3 + len(current_targets))
        
        # 绘制每一层
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
                
                # 添加相关系数文本
                text_angle_rad = theta_rad[j]
                text_radius = r_inner + 0.45
                x = text_radius * np.cos(text_angle_rad)
                y = text_radius * np.sin(text_angle_rad)
                
                val = values.iloc[j]
                if not np.isnan(val):
                    sig_marker = '*' if sig_values.iloc[j] else ''
                    text_val = f'{val:.2f}{sig_marker}'
                    
                    rot = theta_deg[j] - 90 if np.cos(text_angle_rad) > -0.01 else theta_deg[j] + 90
                    
                    ax.text(
                        x, y, text_val,
                        ha='center', va='center',
                        fontsize=4.5, rotation=rot,
                        color='white' if abs(val) > 0.6 else 'black'
                    )
        
        # 组内环形分隔线
        for r in radii:
            circ = Circle((0, 0), r + 0.95, fill=False, edgecolor='#EDEDED', linewidth=0.5, alpha=0.9)
            ax.add_patch(circ)
        
        # 添加特征变量标签（最外圈）
        label_radius = radii.max() + 1.8
        for i in range(len(features)):
            text_angle_rad = theta_rad[i]
            x = label_radius * np.cos(text_angle_rad)
            y = label_radius * np.sin(text_angle_rad)
            rot = theta_deg[i] if np.cos(text_angle_rad) > -0.01 else theta_deg[i] + 180
            
            # 简化变量名
            label_text = features[i].replace('_avg', '').replace('_sum', '').replace('_wsum', '')
            
            ax.text(
                x, y, label_text,
                ha='center', va='center',
                fontsize=3.2, rotation=rot,
                color=current_group_color,
                fontweight='bold'
            )
        
        # 添加分组标签
        group_label_angle_deg = (start_angle_deg + end_angle_deg) / 2
        group_label_angle_rad = np.deg2rad(group_label_angle_deg)
        group_label_radius = radii.max() + 4.7
        x = group_label_radius * np.cos(group_label_angle_rad)
        y = group_label_radius * np.sin(group_label_angle_rad)
        
        ax.text(
            x, y, group_name,
            ha='center', va='center',
            fontsize=18, fontweight='bold',
            color=current_group_color,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                     edgecolor=current_group_color, linewidth=2)
        )

        # 分组径向分隔线（起止角）
        for ang in (start_angle_deg, end_angle_deg):
            a = np.deg2rad(ang)
            x_end = (group_label_radius + 0.5) * np.cos(a)
            y_end = (group_label_radius + 0.5) * np.sin(a)
            ax.plot([0, x_end], [0, y_end], color='#DDDDDD', linewidth=1.0, zorder=0)
    
    # 创建图例（目标分组）
    legend_positions = [
        {'bbox_to_anchor': (0.5, 0.5), 'loc': 'lower right'},
        {'bbox_to_anchor': (0.5, 0.5), 'loc': 'lower left'},
        {'bbox_to_anchor': (0.5, 0.5), 'loc': 'upper right'},
        {'bbox_to_anchor': (0.5, 0.5), 'loc': 'upper left'}
    ]
    
    for i, group_name in enumerate(group_names_list):
        handles_for_group = []
        current_targets = all_target_names[group_name]
        current_group_color = group_legend_colors[i]
        
        for j, target_name in enumerate(current_targets):
            marker_shapes = ['o', 's', '^', 'v', 'D', '*', 'p', 'h']
            marker = marker_shapes[j % len(marker_shapes)]
            
            handle = Line2D(
                [0], [0], marker=marker, color='w',
                markerfacecolor=current_group_color,
                markersize=10, label=target_name
            )
            handles_for_group.append(handle)
        
        if i < len(legend_positions):
            leg = ax.legend(
                handles=handles_for_group,
                title=f"{group_name} correlates with:",
                **legend_positions[i],
                fontsize=10,
                title_fontsize=12,
                frameon=True,
                fancybox=True,
                shadow=True
            )
            leg.get_title().set_fontweight('bold')
            ax.add_artist(leg)
    
    # 添加颜色条
    cax = fig.add_axes([0.2, 0.08, 0.6, 0.012])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(f"{selected_method.capitalize()} correlation", size=16, weight='bold')
    cbar.ax.tick_params(size=12, labelsize=12)
    cbar.set_ticks([-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0])
    
    # 设置坐标轴范围
    max_radius = max([max(np.arange(3, 3 + len(all_target_names[g]))) 
                     for g in group_names_list]) + 6
    ax.set_xlim(-max_radius, max_radius)
    ax.set_ylim(-max_radius, max_radius)
    
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
    start_angle = i * angle_per_group + gap/2
    end_angle = (i + 1) * angle_per_group - gap/2
    marker_angle = (start_angle + end_angle) / 2
    
    sector_params[group_name] = {
        'start': start_angle,
        'end': end_angle,
        'marker_angle': marker_angle
    }

# 准备目标字典（每组对应其他组）
all_target_names = {}
for group_name in group_names:
    all_target_names[group_name] = [g for g in group_names if g != group_name]

# 绘制
fig = create_full_ring_plot(
    all_data=all_correlation_data,
    all_feature_names=var_groups,
    all_target_names=all_target_names,
    color_palette=select_color,
    sector_params=sector_params
)

# =============================================================================
# 9. 保存结果
# =============================================================================
print("\n步骤 5/6: 保存图表...")

# 保存PNG（高分辨率）
png_path = os.path.join(output_directory, "petal_correlation_plot.png")
fig.savefig(png_path, dpi=1200, bbox_inches='tight', facecolor='white')
print(f"  - 已保存: {png_path}")

# 保存PDF（矢量图）
pdf_path = os.path.join(output_directory, "petal_correlation_plot.pdf")
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"  - 已保存: {pdf_path}")

# 额外导出：将每个变量组的相关性矩阵与显著性矩阵保存为CSV
csv_dir = os.path.join(data_directory, "petal_tables")
os.makedirs(csv_dir, exist_ok=True)
for group_name in group_names:
    corr_df = all_correlation_data[group_name]['correlation_df']
    sig_df = all_correlation_data[group_name]['p_value_df']
    corr_csv = os.path.join(csv_dir, f"{group_name}_correlation.csv")
    sig_csv = os.path.join(csv_dir, f"{group_name}_significance.csv")
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
print(f"  - PDF: {pdf_path}")
print(f"  - CSV表: {csv_dir}")

print("\n版式: Arial 字体, 蓝-灰-黄配色, 1200 dpi, 单图导出")
