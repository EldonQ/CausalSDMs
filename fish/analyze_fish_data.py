# -*- coding: utf-8 -*-
"""
DASCO_AlienCoordinates_SInAS_2.4.1.csv 数据分析脚本
===================================================
分析入侵鱼类分布数据，包含物种分类、地理分布、生物群系、数据库来源等维度。

数据来源: DASCO项目，OBIS等数据库
文件大小: ~2.24GB
列字段: taxon, location, Realm, Longitude, Latitude, Database, occurrenceID, eventDate

Author: AI Analysis
Date: 2026-03-25
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 非GUI后端
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# ============ 配置参数 ============
CSV_PATH = r'E:\CausalSDMs\fish\DASCO_AlienCoordinates_SInAS_2.4.1.csv'
OUTPUT_DIR = r'E:\CausalSDMs\fish\analysis_output'
CHUNK_SIZE = 500_000  # 每块读取50万行

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============ 1. 数据概览 ============
def step1_data_overview():
    """步骤1: 数据基本信息概览"""
    print("\n" + "=" * 70)
    print("【步骤1】数据基本信息概览")
    print("=" * 70)

    # 先读取少量数据获取列信息
    df_sample = pd.read_csv(CSV_PATH, nrows=1000, encoding='utf-8')
    print(f"\n列名: {list(df_sample.columns)}")
    print(f"列数: {len(df_sample.columns)}")

    # 数据类型
    print("\n数据类型:")
    print(df_sample.dtypes)

    # 估算总行数（通过文件大小估算）
    file_size_gb = os.path.getsize(CSV_PATH) / (1024 ** 3)
    avg_row_size_kb = df_sample.memory_usage(deep=True).sum() / len(df_sample) / 1024
    estimated_rows = int(file_size_gb * 1024 * 1024 / avg_row_size_kb)
    print(f"\n文件大小: {file_size_gb:.2f} GB")
    print(f"估算总行数: 约 {estimated_rows:,} 行")

    # 空值示例
    print("\n前5行数据预览:")
    print(df_sample.head())


# ============ 2. 分块读取统计 ============
def step2_chunked_statistics():
    """步骤2: 分块读取并统计所有数值"""
    print("\n" + "=" * 70)
    print("【步骤2】分块读取并统计")
    print("=" * 70)

    # 初始化统计变量
    total_rows = 0
    taxon_counts = {}
    location_counts = {}
    realm_counts = {}
    database_counts = {}
    longitudes = []
    latitudes = []
    event_dates = []
    null_counts = {col: 0 for col in ['taxon', 'location', 'Realm', 'Longitude', 'Latitude', 'Database', 'occurrenceID', 'eventDate']}

    # 分类统计
    realms_in_data = {}  # {realm: {taxon_count, record_count}}

    print("开始分块读取数据...")
    chunk_num = 0

    for chunk in pd.read_csv(CSV_PATH, encoding='utf-8', chunksize=CHUNK_SIZE,
                             dtype={'Longitude': 'float64', 'Latitude': 'float64'}):
        chunk_num += 1
        n = len(chunk)
        total_rows += n

        # 空值统计
        for col in null_counts:
            null_counts[col] += chunk[col].isna().sum()

        # 物种统计
        for taxon, cnt in chunk['taxon'].value_counts().items():
            taxon_counts[taxon] = taxon_counts.get(taxon, 0) + cnt

        # 地理位置统计
        for loc, cnt in chunk['location'].value_counts().items():
            location_counts[loc] = location_counts.get(loc, 0) + cnt

        # 生物群系统计
        for realm, cnt in chunk['Realm'].value_counts().items():
            realm_counts[realm] = realm_counts.get(realm, 0) + cnt
            # 记录每个realm下的物种
            if realm not in realms_in_data:
                realms_in_data[realm] = {'taxons': set(), 'records': 0}
            realms_in_data[realm]['taxons'].update(chunk['taxon'].dropna().unique())
            realms_in_data[realm]['records'] += cnt

        # 数据库来源统计
        for db, cnt in chunk['Database'].value_counts().items():
            database_counts[db] = database_counts.get(db, 0) + cnt

        # 经纬度
        longitudes.extend(chunk['Longitude'].dropna().tolist())
        latitudes.extend(chunk['Latitude'].dropna().tolist())

        # 事件日期
        event_dates.extend(chunk['eventDate'].dropna().tolist())

        print(f"  已处理 {chunk_num} 块 ({total_rows:,} 行)...")

    print(f"\n总记录数: {total_rows:,}")

    # 保存统计结果
    results = {
        'total_rows': total_rows,
        'taxon_counts': taxon_counts,
        'location_counts': location_counts,
        'realm_counts': realm_counts,
        'database_counts': database_counts,
        'longitudes': longitudes,
        'latitudes': latitudes,
        'event_dates': event_dates,
        'null_counts': null_counts,
        'realms_in_data': realms_in_data,
    }

    return results


# ============ 3. 生成分析报告 ============
def step3_analysis_report(results):
    """步骤3: 生成详细分析报告"""
    print("\n" + "=" * 70)
    print("【步骤3】详细分析报告")
    print("=" * 70)

    total = results['total_rows']

    # 3.1 物种分析
    print("\n--- 3.1 物种(Taxon)分析 ---")
    print(f"物种总数: {len(results['taxon_counts'])}")
    print("\n记录数Top20物种:")
    sorted_taxons = sorted(results['taxon_counts'].items(), key=lambda x: -x[1])
    for i, (taxon, cnt) in enumerate(sorted_taxons[:20], 1):
        print(f"  {i:2d}. {taxon}: {cnt:,} 条记录 ({cnt/total*100:.2f}%)")

    # 3.2 地理位置分析
    print("\n--- 3.2 地理位置(Location)分析 ---")
    print(f"不同地理位置数量: {len(results['location_counts'])}")
    print("\n记录数Top20地理位置:")
    sorted_locations = sorted(results['location_counts'].items(), key=lambda x: -x[1])
    for i, (loc, cnt) in enumerate(sorted_locations[:20], 1):
        print(f"  {i:2d}. {loc}: {cnt:,} 条记录 ({cnt/total*100:.2f}%)")

    # 3.3 生物群系分析
    print("\n--- 3.3 生物群系(Realm)分析 ---")
    for realm, data in sorted(results['realm_counts'].items(), key=lambda x: -x[1]):
        taxons = results['realms_in_data'].get(realm, {}).get('taxons', set())
        print(f"  {realm}: {data:,} 条记录, {len(taxons)} 个物种 ({data/total*100:.2f}%)")

    # 3.4 数据库来源分析
    print("\n--- 3.4 数据库(Database)来源分析 ---")
    for db, cnt in sorted(results['database_counts'].items(), key=lambda x: -x[1]):
        print(f"  {db}: {cnt:,} 条记录 ({cnt/total*100:.2f}%)")

    # 3.5 经纬度范围
    print("\n--- 3.5 地理坐标范围 ---")
    lons = np.array(results['longitudes'])
    lats = np.array(results['latitudes'])
    print(f"  经度范围: [{lons.min():.4f}, {lons.max():.4f}]")
    print(f"  纬度范围: [{lats.min():.4f}, {lats.max():.4f}]")
    print(f"  经度均值: {lons.mean():.4f}, 标准差: {lons.std():.4f}")
    print(f"  纬度均值: {lats.mean():.4f}, 标准差: {lats.std():.4f}")
    print(f"  有效坐标数: {len(lons):,}")

    # 3.6 空值统计
    print("\n--- 3.6 空值(Null)统计 ---")
    for col, cnt in results['null_counts'].items():
        print(f"  {col}: {cnt:,} 个空值 ({cnt/total*100:.2f}%)")

    # 3.7 时间分析
    print("\n--- 3.7 事件日期(eventDate)分析 ---")
    if results['event_dates']:
        valid_dates = [d for d in results['event_dates'] if str(d).strip() != '']
        print(f"  有效日期记录数: {len(valid_dates):,} ({len(valid_dates)/total*100:.2f}%)")
        if valid_dates:
            # 尝试解析日期
            try:
                dates_series = pd.to_datetime(valid_dates, errors='coerce')
                valid_dt = dates_series.dropna()
                if len(valid_dt) > 0:
                    print(f"  最早日期: {valid_dt.min()}")
                    print(f"  最晚日期: {valid_dt.max()}")
                    print(f"  日期范围跨度: {(valid_dt.max() - valid_dt.min()).days} 天")
            except:
                print(f"  日期样本: {valid_dates[:5]}")
    else:
        print("  无有效日期数据")

    return sorted_taxons, sorted_locations


# ============ 4. 生成可视化 ============
def step4_generate_visualizations(results):
    """步骤4: 生成可视化图表"""
    print("\n" + "=" * 70)
    print("【步骤4】生成可视化图表")
    print("=" * 70)

    total = results['total_rows']
    sorted_taxons = sorted(results['taxon_counts'].items(), key=lambda x: -x[1])
    sorted_locations = sorted(results['location_counts'].items(), key=lambda x: -x[1])

    # 图1: 物种分布 - Top15柱状图
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # 1.1 物种Top15
    ax1 = axes[0, 0]
    top15_taxons = sorted_taxons[:15]
    names = [t[0][:30] for t in top15_taxons]
    values = [t[1] for t in top15_taxons]
    bars = ax1.barh(range(len(names)), values, color='steelblue')
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=8)
    ax1.invert_yaxis()
    ax1.set_xlabel('Record Count')
    ax1.set_title(f'Top 15 Species by Record Count\n(Total: {total:,} records, {len(results["taxon_counts"])} species)', fontsize=10)
    for i, v in enumerate(values):
        ax1.text(v + max(values)*0.01, i, f'{v:,}', va='center', fontsize=7)

    # 1.2 生物群系饼图
    ax2 = axes[0, 1]
    realm_data = sorted(results['realm_counts'].items(), key=lambda x: -x[1])
    labels = [r[0] for r in realm_data]
    sizes = [r[1] for r in realm_data]
    colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
    wedges, texts, autotexts = ax2.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
    ax2.set_title(f'Distribution by Realm\n({len(labels)} realms)', fontsize=10)

    # 1.3 数据库来源柱状图
    ax3 = axes[1, 0]
    db_data = sorted(results['database_counts'].items(), key=lambda x: -x[1])
    db_names = [d[0] for d in db_data]
    db_values = [d[1] for d in db_data]
    bars3 = ax3.bar(range(len(db_names)), db_values, color='coral')
    ax3.set_xticks(range(len(db_names)))
    ax3.set_xticklabels(db_names, rotation=30, ha='right', fontsize=9)
    ax3.set_ylabel('Record Count')
    ax3.set_title(f'Distribution by Database\n({len(db_names)} sources)', fontsize=10)
    for i, v in enumerate(db_values):
        ax3.text(i, v + max(db_values)*0.01, f'{v:,}', ha='center', fontsize=8)

    # 1.4 地理坐标散点图
    ax4 = axes[1, 1]
    lons = np.array(results['longitudes'])
    lats = np.array(results['latitudes'])
    # 随机采样绘制（如果数据量太大）
    sample_size = min(100000, len(lons))
    if sample_size < len(lons):
        indices = np.random.choice(len(lons), sample_size, replace=False)
        plot_lons, plot_lats = lons[indices], lats[indices]
        title_suffix = f' (sampled {sample_size:,} of {len(lons):,})'
    else:
        plot_lons, plot_lats = lons, lats
        title_suffix = f' ({len(lons):,} points)'

    scatter = ax4.scatter(plot_lons, plot_lats, c='blue', alpha=0.1, s=1)
    ax4.set_xlabel('Longitude')
    ax4.set_ylabel('Latitude')
    ax4.set_title(f'Geographic Distribution of Records{title_suffix}', fontsize=10)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    fig1_path = os.path.join(OUTPUT_DIR, 'fig1_overview.png')
    plt.savefig(fig1_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: {fig1_path}")

    # 图2: 地理位置Top20 + 经纬度分布
    fig2, axes2 = plt.subplots(1, 2, figsize=(16, 8))

    # 2.1 地理位置Top20
    ax5 = axes2[0]
    top20_locs = sorted_locations[:20]
    loc_names = [l[0][:35] for l in top20_locs]
    loc_values = [l[1] for l in top20_locs]
    bars5 = ax5.barh(range(len(loc_names)), loc_values, color='seagreen')
    ax5.set_yticks(range(len(loc_names)))
    ax5.set_yticklabels(loc_names, fontsize=8)
    ax5.invert_yaxis()
    ax5.set_xlabel('Record Count')
    ax5.set_title(f'Top 20 Locations by Record Count\n({len(results["location_counts"])} unique locations)', fontsize=10)
    for i, v in enumerate(loc_values):
        ax5.text(v + max(loc_values)*0.01, i, f'{v:,}', va='center', fontsize=7)

    # 2.2 经纬度直方图
    ax6 = axes2[1]
    ax6.hist(lons, bins=50, alpha=0.6, label='Longitude', color='blue')
    ax6.hist(lats, bins=50, alpha=0.6, label='Latitude', color='orange')
    ax6.set_xlabel('Coordinate Value')
    ax6.set_ylabel('Frequency')
    ax6.set_title(f'Distribution of Coordinates\n(Lon range: [{lons.min():.2f}, {lons.max():.2f}], Lat range: [{lats.min():.2f}, {lats.max():.2f}])', fontsize=10)
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    plt.tight_layout()
    fig2_path = os.path.join(OUTPUT_DIR, 'fig2_locations_coordinates.png')
    plt.savefig(fig2_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: {fig2_path}")

    # 图3: 各生物群系物种数量和记录数
    fig3, ax7 = plt.subplots(figsize=(12, 6))
    realms_sorted = sorted(results['realms_in_data'].items(),
                          key=lambda x: -results['realm_counts'].get(x[0], 0))
    realm_names = [r[0] for r in realms_sorted]
    taxon_counts = [len(r[1]['taxons']) for r in realms_sorted]
    record_counts = [results['realm_counts'].get(r[0], 0) for r in realms_sorted]

    x = np.arange(len(realm_names))
    width = 0.35
    bars1 = ax7.bar(x - width/2, taxon_counts, width, label='Species Count', color='steelblue')
    ax7_twin = ax7.twinx()
    bars2 = ax7_twin.bar(x + width/2, record_counts, width, label='Record Count', color='coral')

    ax7.set_xticks(x)
    ax7.set_xticklabels(realm_names, rotation=30, ha='right', fontsize=9)
    ax7.set_ylabel('Species Count', color='steelblue')
    ax7_twin.set_ylabel('Record Count', color='coral')
    ax7.set_title('Species and Record Count by Realm', fontsize=12)
    ax7.legend(loc='upper left')
    ax7_twin.legend(loc='upper right')

    plt.tight_layout()
    fig3_path = os.path.join(OUTPUT_DIR, 'fig3_realm_species_records.png')
    plt.savefig(fig3_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: {fig3_path}")

    # 图4: 热力图 - 按经纬度网格统计
    fig4, ax8 = plt.subplots(figsize=(14, 8))
    if len(lons) > 1000:
        # 创建2D直方图作为热力图
        lon_bins = np.linspace(lons.min(), lons.max(), 100)
        lat_bins = np.linspace(lats.min(), lats.max(), 100)
        h, xedges, yedges = np.histogram2d(lons, lats, bins=[lon_bins, lat_bins])
        h = np.log1p(h)  # 对数变换便于可视化
        im = ax8.imshow(h.T, origin='lower',
                        extent=[lons.min(), lons.max(), lats.min(), lats.max()],
                        aspect='auto', cmap='YlOrRd')
        plt.colorbar(im, ax=ax8, label='Log(Count+1)')
    ax8.set_xlabel('Longitude')
    ax8.set_ylabel('Latitude')
    ax8.set_title(f'Spatial Density Heatmap (log scale)\n({len(lons):,} records)', fontsize=12)
    ax8.grid(True, alpha=0.2)

    plt.tight_layout()
    fig4_path = os.path.join(OUTPUT_DIR, 'fig4_spatial_heatmap.png')
    plt.savefig(fig4_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: {fig4_path}")

    print("\n所有图表已生成完毕！")


# ============ 5. 导出详细数据 ============
def step5_export_data(results):
    """步骤5: 导出详细统计数据到CSV"""
    print("\n" + "=" * 70)
    print("【步骤5】导出详细统计数据")
    print("=" * 70)

    total = results['total_rows']

    # 5.1 物种统计
    df_taxons = pd.DataFrame([
        {'taxon': t[0], 'record_count': t[1], 'percentage': t[1]/total*100}
        for t in sorted(results['taxon_counts'].items(), key=lambda x: -x[1])
    ])
    taxon_csv = os.path.join(OUTPUT_DIR, 'stats_taxons.csv')
    df_taxons.to_csv(taxon_csv, index=False, encoding='utf-8-sig')
    print(f"  已导出: {taxon_csv} ({len(df_taxons)} 条)")

    # 5.2 地理位置统计
    df_locations = pd.DataFrame([
        {'location': l[0], 'record_count': l[1], 'percentage': l[1]/total*100}
        for l in sorted(results['location_counts'].items(), key=lambda x: -x[1])
    ])
    loc_csv = os.path.join(OUTPUT_DIR, 'stats_locations.csv')
    df_locations.to_csv(loc_csv, index=False, encoding='utf-8-sig')
    print(f"  已导出: {loc_csv} ({len(df_locations)} 条)")

    # 5.3 生物群系统计
    realm_rows = []
    for realm, cnt in results['realm_counts'].items():
        taxons = results['realms_in_data'].get(realm, {}).get('taxons', set())
        realm_rows.append({
            'realm': realm,
            'record_count': cnt,
            'species_count': len(taxons),
            'percentage': cnt/total*100
        })
    df_realms = pd.DataFrame(realm_rows).sort_values('record_count', ascending=False)
    realm_csv = os.path.join(OUTPUT_DIR, 'stats_realms.csv')
    df_realms.to_csv(realm_csv, index=False, encoding='utf-8-sig')
    print(f"  已导出: {realm_csv} ({len(df_realms)} 条)")

    # 5.4 数据库统计
    df_databases = pd.DataFrame([
        {'database': d[0], 'record_count': d[1], 'percentage': d[1]/total*100}
        for d in sorted(results['database_counts'].items(), key=lambda x: -x[1])
    ])
    db_csv = os.path.join(OUTPUT_DIR, 'stats_databases.csv')
    df_databases.to_csv(db_csv, index=False, encoding='utf-8-sig')
    print(f"  已导出: {db_csv} ({len(df_databases)} 条)")

    # 5.5 综合摘要报告
    summary_path = os.path.join(OUTPUT_DIR, 'analysis_summary.txt')
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("DASCO 入侵鱼类数据 综合分析报告\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"数据文件: {CSV_PATH}\n")
        f.write(f"总记录数: {total:,}\n")
        f.write(f"物种数量: {len(results['taxon_counts'])}\n")
        f.write(f"地理位置数量: {len(results['location_counts'])}\n")
        f.write(f"生物群系数量: {len(results['realm_counts'])}\n")
        f.write(f"数据库来源数量: {len(results['database_counts'])}\n\n")

        lons = np.array(results['longitudes'])
        lats = np.array(results['latitudes'])
        f.write("地理范围:\n")
        f.write(f"  经度: [{lons.min():.4f}, {lons.max():.4f}]\n")
        f.write(f"  纬度: [{lats.min():.4f}, {lats.max():.4f}]\n\n")

        f.write("空值统计:\n")
        for col, cnt in results['null_counts'].items():
            f.write(f"  {col}: {cnt:,} ({cnt/total*100:.2f}%)\n")

        f.write("\n各生物群系详情:\n")
        for realm, data in sorted(results['realm_counts'].items(), key=lambda x: -x[1]):
            taxons = results['realms_in_data'].get(realm, {}).get('taxons', set())
            f.write(f"  {realm}: {data:,} 条记录, {len(taxons)} 个物种\n")

        f.write("\nTop10 物种:\n")
        for i, (taxon, cnt) in enumerate(sorted(results['taxon_counts'].items(), key=lambda x: -x[1])[:10], 1):
            f.write(f"  {i:2d}. {taxon}: {cnt:,} 条 ({cnt/total*100:.2f}%)\n")

    print(f"  已导出: {summary_path}")
    print("\n所有数据导出完毕！")


# ============ 主函数 ============
def main():
    print("=" * 70)
    print("DASCO 入侵鱼类数据 (DASCO_AlienCoordinates_SInAS_2.4.1.csv) 分析脚本")
    print("=" * 70)

    # 步骤1: 数据概览
    step1_data_overview()

    # 步骤2: 分块统计（耗时最长）
    results = step2_chunked_statistics()

    # 步骤3: 分析报告
    step3_analysis_report(results)

    # 步骤4: 可视化
    step4_generate_visualizations(results)

    # 步骤5: 导出数据
    step5_export_data(results)

    print("\n" + "=" * 70)
    print("分析完成！所有结果保存在:")
    print(f"  {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
