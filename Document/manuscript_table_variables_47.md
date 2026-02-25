### 论文用环境变量清单（47项）

说明：本表基于 `scripts/variables_selected_47.csv` 列表生成，单位与变量含义参考 EarthEnv-Streams 原始说明文档（ReadMe）。若原始栅格采用缩放存储（如温度以 °C×10、坡度以 度×100、土壤 pH 以 ×10），则在“单位（原始）”中明确标注。

| 变量名 | 文件 | 波段 | 分组 | 单位（原始） | 中文描述 |
|---|---|---:|---|---|---|
| dem_range | elevation.tif | 3 | G1_TopoSlopeFlow | m | 高程范围（最大值 − 最小值） |
| dem_avg | elevation.tif | 4 | G1_TopoSlopeFlow | m | 平均高程 |
| slope_range | slope.tif | 3 | G1_TopoSlopeFlow | 度×100 | 坡度范围（最大值 − 最小值），以度×100存储 |
| slope_avg | slope.tif | 4 | G1_TopoSlopeFlow | 度×100 | 平均坡度，以度×100存储 |
| flow_length | flow_acc.tif | 1 | G1_TopoSlopeFlow | 计数 | 上游河道栅格单元数量（count） |
| flow_acc | flow_acc.tif | 2 | G1_TopoSlopeFlow | 计数 | 上游汇水分区栅格单元数量（count） |
| hydro_wavg_01 | hydroclim_weighted_average+sum.tif | 1 | G2_Hydroclim_wavg | °C×10 | 年均上游温度（Bioclim 1） |
| hydro_wavg_02 | hydroclim_weighted_average+sum.tif | 2 | G2_Hydroclim_wavg | °C×10 | 上游月昼夜温差均值（各月 Tmax−Tmin 的月均；Bioclim 2） |
| hydro_wavg_03 | hydroclim_weighted_average+sum.tif | 3 | G2_Hydroclim_wavg | ×100 | 上游等温性（Bioclim 3 = 02/07 ×100；无量纲） |
| hydro_wavg_04 | hydroclim_weighted_average+sum.tif | 4 | G2_Hydroclim_wavg | °C×10 | 上游温度季节性（标准差×100，Bioclim 4；原始温度以 °C×10 存储） |
| hydro_wavg_05 | hydroclim_weighted_average+sum.tif | 5 | G2_Hydroclim_wavg | °C×10 | 最暖月上游最高温（Bioclim 5） |
| hydro_wavg_06 | hydroclim_weighted_average+sum.tif | 6 | G2_Hydroclim_wavg | °C×10 | 最冷月上游最低温（Bioclim 6） |
| hydro_wavg_07 | hydroclim_weighted_average+sum.tif | 7 | G2_Hydroclim_wavg | °C×10 | 上游年温差（Bioclim 7 = 05−06） |
| hydro_wavg_08 | hydroclim_weighted_average+sum.tif | 8 | G2_Hydroclim_wavg | °C×10 | 最湿季上游平均温（Bioclim 8） |
| hydro_wavg_09 | hydroclim_weighted_average+sum.tif | 9 | G2_Hydroclim_wavg | °C×10 | 最干季上游平均温（Bioclim 9） |
| hydro_wavg_10 | hydroclim_weighted_average+sum.tif | 10 | G2_Hydroclim_wavg | °C×10 | 最暖季上游平均温（Bioclim 10） |
| hydro_wavg_11 | hydroclim_weighted_average+sum.tif | 11 | G2_Hydroclim_wavg | °C×10 | 最冷季上游平均温（Bioclim 11） |
| hydro_wavg_12 | hydroclim_weighted_average+sum.tif | 12 | G2_Hydroclim_wavg | mm | 年上游降水量（Bioclim 12） |
| hydro_wavg_13 | hydroclim_weighted_average+sum.tif | 13 | G2_Hydroclim_wavg | mm | 最湿月上游降水量（Bioclim 13） |
| hydro_wavg_14 | hydroclim_weighted_average+sum.tif | 14 | G2_Hydroclim_wavg | mm | 最干月上游降水量（Bioclim 14） |
| hydro_wavg_15 | hydroclim_weighted_average+sum.tif | 15 | G2_Hydroclim_wavg | ×100 | 上游降水季节性（Bioclim 15，变异系数×100；无量纲） |
| hydro_wavg_16 | hydroclim_weighted_average+sum.tif | 16 | G2_Hydroclim_wavg | mm | 最湿季上游降水量（Bioclim 16） |
| hydro_wavg_17 | hydroclim_weighted_average+sum.tif | 17 | G2_Hydroclim_wavg | mm | 最干季上游降水量（Bioclim 17） |
| hydro_wavg_18 | hydroclim_weighted_average+sum.tif | 18 | G2_Hydroclim_wavg | mm | 最暖季上游降水量（Bioclim 18） |
| hydro_wavg_19 | hydroclim_weighted_average+sum.tif | 19 | G2_Hydroclim_wavg | mm | 最冷季上游降水量（Bioclim 19） |
| lc_wavg_01 | landcover_weighted_average.tif | 1 | G3_Landcover_wavg | % | 上游土地覆盖：常绿/落叶针叶林（Evergreen/deciduous needleleaf trees）加权平均占比 |
| lc_wavg_02 | landcover_weighted_average.tif | 2 | G3_Landcover_wavg | % | 上游土地覆盖：常绿阔叶林（Evergreen broadleaf trees）加权平均占比 |
| lc_wavg_03 | landcover_weighted_average.tif | 3 | G3_Landcover_wavg | % | 上游土地覆盖：落叶阔叶林（Deciduous broadleaf trees）加权平均占比 |
| lc_wavg_04 | landcover_weighted_average.tif | 4 | G3_Landcover_wavg | % | 上游土地覆盖：混交/其他乔木（Mixed/other trees）加权平均占比 |
| lc_wavg_05 | landcover_weighted_average.tif | 5 | G3_Landcover_wavg | % | 上游土地覆盖：灌丛（Shrubs）加权平均占比 |
| lc_wavg_06 | landcover_weighted_average.tif | 6 | G3_Landcover_wavg | % | 上游土地覆盖：草本植被（Herbaceous vegetation）加权平均占比 |
| lc_wavg_07 | landcover_weighted_average.tif | 7 | G3_Landcover_wavg | % | 上游土地覆盖：耕地/管理植被（Cultivated and managed vegetation）加权平均占比 |
| lc_wavg_08 | landcover_weighted_average.tif | 8 | G3_Landcover_wavg | % | 上游土地覆盖：常年淹水灌丛/草本（Regularly flooded shrub/herbaceous vegetation）加权平均占比 |
| lc_wavg_09 | landcover_weighted_average.tif | 9 | G3_Landcover_wavg | % | 上游土地覆盖：城市/建成区（Urban/built-up）加权平均占比 |
| lc_wavg_10 | landcover_weighted_average.tif | 10 | G3_Landcover_wavg | % | 上游土地覆盖：雪/冰（Snow/ice）加权平均占比 |
| lc_wavg_11 | landcover_weighted_average.tif | 11 | G3_Landcover_wavg | % | 上游土地覆盖：裸地/稀疏植被（Barren lands/sparse vegetation）加权平均占比 |
| lc_wavg_12 | landcover_weighted_average.tif | 12 | G3_Landcover_wavg | % | 上游土地覆盖：开阔水体（Open water）加权平均占比 |
| soil_wavg_01 | soil_weighted_average.tif | 1 | G4_Soil_wavg | g/kg | 上游土壤有机碳（Soil organic carbon）加权平均 |
| soil_wavg_02 | soil_weighted_average.tif | 2 | G4_Soil_wavg | pH×10 | 上游土壤 pH（H2O）加权平均（以 pH×10 存储） |
| soil_wavg_03 | soil_weighted_average.tif | 3 | G4_Soil_wavg | % | 上游土壤砂含量（Sand content）加权平均 |
| soil_wavg_04 | soil_weighted_average.tif | 4 | G4_Soil_wavg | % | 上游土壤粉砂含量（Silt content）加权平均 |
| soil_wavg_05 | soil_weighted_average.tif | 5 | G4_Soil_wavg | % | 上游土壤黏土含量（Clay content）加权平均 |
| soil_wavg_06 | soil_weighted_average.tif | 6 | G4_Soil_wavg | % | 上游粗颗粒体积分数（>2 mm；Coarse fragments）加权平均 |
| soil_wavg_07 | soil_weighted_average.tif | 7 | G4_Soil_wavg | cmol/kg | 上游阳离子交换量（Cation exchange capacity）加权平均 |
| soil_wavg_08 | soil_weighted_average.tif | 8 | G4_Soil_wavg | kg/m³ | 上游细土体密度（Bulk density of fine earth）加权平均 |
| soil_wavg_09 | soil_weighted_average.tif | 9 | G4_Soil_wavg | cm | 基岩埋深（R层；至最大240 cm）上游加权平均 |
| soil_wavg_10 | soil_weighted_average.tif | 10 | G4_Soil_wavg | % | R层出现概率（0–100%）上游加权平均 |

注：依据项目脚本约定，分析环节对若干变量做单位还原以便解释与制图：温度类（hydro_wavg_01–11）÷10 转为 °C；坡度（slope_*）÷100 转为度；土壤 pH（soil_wavg_02）÷10 转为标准 pH 值。


