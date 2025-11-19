## EarthEnv-Streams 1 km 环境变量（324 层）中文详解

本文件系统汇总 `EarthEnv-Streams` 数据集中 324 个近全球（30 弧秒 ≈ 1 km）淡水相关环境变量，逐一给出变量代码、英文含义、中文解释、聚合方式与单位，方便科研引用与建模选取。

- Extent: (-145, -56, 180, 60)
- Cell size: 0.0083333333° (30 arcsec)
- Projection: WGS84

引用与来源（务必在论文中引用）：
- Domisch, S., Amatulli, G., Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution. Scientific Data.
- HydroSHEDS（Lehner et al., 2008）、GLWD（Lehner & Döll, 2004）、WorldClim（Hijmans et al., 2005）、CLC（Tuanmu & Jetz, 2014）、美国地质调查局地表地质 USGS（见 ReadMe 原文链接）、ISRIC SoilGrids（Hengl et al., 2014）。

说明：
- 坡度单位“度×100”表示原始坡度（°）乘以 100 的整型存储。
- 气温单位“℃×10”表示原始气温（℃）乘以 10 的整型存储。
- “加权”表示按照上游汇流贡献在河网中进行距离与权重传播后的加权统计。

---

### 目录（按数据文件分组）
- 地形：`elevation.nc`（4）
- 坡度：`slope.nc`（4）
- 汇流累积：`flow_acc.nc`（2）
- 月尺度气候（上游平均/和）：`monthly_tmin_average.nc`（12），`monthly_tmax_average.nc`（12），`monthly_prec_sum.nc`（12）
- 月尺度气候（上游加权平均/和）：`monthly_tmin_weighted_average.nc`（12），`monthly_tmax_weighted_average.nc`（12），`monthly_prec_weighted_sum.nc`（12）
- 上游 Bioclim（平均/和）：`hydroclim_average+sum.nc`（19）
- 上游 Bioclim（加权平均/和）：`hydroclim_weighted_average+sum.nc`（19）
- 土地覆盖（最小/最大/范围/平均/加权平均）：`landcover_*.nc`（共 60）
- 地表地质（加权计数）：`geology_weighted_sum.nc`（92）
- 土壤（最小/最大/范围/平均/加权平均）：`soil_*.nc`（共 50）
- 质量控制：`quality_control.nc`（2）

---

### 1) 地形 Elevation — `elevation.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | dem_min | Minimum elevation | 最小高程 | minimum | m |
| 2 | dem_max | Maximum elevation | 最大高程 | maximum | m |
| 3 | dem_range | Elevation range | 高程范围（最大-最小） | range | m |
| 4 | dem_avg | Average elevation | 平均高程 | average | m |

来源：HydroSHEDS。

---

### 2) 坡度 Slope — `slope.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | slope_min | Minimum slope | 最小坡度 | minimum | 度×100 |
| 2 | slope_max | Maximum slope | 最大坡度 | maximum | 度×100 |
| 3 | slope_range | Slope range | 坡度范围（最大-最小） | range | 度×100 |
| 4 | slope_avg | Average slope | 平均坡度 | average | 度×100 |

来源：HydroSHEDS。

---

### 3) 汇流累积 Flow Accumulation — `flow_acc.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | flow_length | Number of upstream stream grid cells | 上游河道网格数量（沿河道长度累计，计数） | sum | count |
| 2 | flow_acc | Number of upstream catchment grid cells | 上游汇水区网格数量（面积累计，计数） | sum | count |

来源：HydroSHEDS。

---

### 4) 月均最低气温（上游平均）— `monthly_tmin_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | tmin_avg_01 | Min. monthly air temperature for January | 1 月上游最低月气温（平均） | average | ℃×10 |
| 2 | tmin_avg_02 | ... February | 2 月上游最低月气温（平均） | average | ℃×10 |
| 3 | tmin_avg_03 | ... March | 3 月上游最低月气温（平均） | average | ℃×10 |
| 4 | tmin_avg_04 | ... April | 4 月上游最低月气温（平均） | average | ℃×10 |
| 5 | tmin_avg_05 | ... May | 5 月上游最低月气温（平均） | average | ℃×10 |
| 6 | tmin_avg_06 | ... June | 6 月上游最低月气温（平均） | average | ℃×10 |
| 7 | tmin_avg_07 | ... July | 7 月上游最低月气温（平均） | average | ℃×10 |
| 8 | tmin_avg_08 | ... August | 8 月上游最低月气温（平均） | average | ℃×10 |
| 9 | tmin_avg_09 | ... September | 9 月上游最低月气温（平均） | average | ℃×10 |
| 10 | tmin_avg_10 | ... October | 10 月上游最低月气温（平均） | average | ℃×10 |
| 11 | tmin_avg_11 | ... November | 11 月上游最低月气温（平均） | average | ℃×10 |
| 12 | tmin_avg_12 | ... December | 12 月上游最低月气温（平均） | average | ℃×10 |

来源：WorldClim。

---

### 5) 月均最高气温（上游平均）— `monthly_tmax_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | tmax_avg_01 | Max. monthly air temperature for January | 1 月上游最高月气温（平均） | average | ℃×10 |
| 2 | tmax_avg_02 | ... February | 2 月上游最高月气温（平均） | average | ℃×10 |
| 3 | tmax_avg_03 | ... March | 3 月上游最高月气温（平均） | average | ℃×10 |
| 4 | tmax_avg_04 | ... April | 4 月上游最高月气温（平均） | average | ℃×10 |
| 5 | tmax_avg_05 | ... May | 5 月上游最高月气温（平均） | average | ℃×10 |
| 6 | tmax_avg_06 | ... June | 6 月上游最高月气温（平均） | average | ℃×10 |
| 7 | tmax_avg_07 | ... July | 7 月上游最高月气温（平均） | average | ℃×10 |
| 8 | tmax_avg_08 | ... August | 8 月上游最高月气温（平均） | average | ℃×10 |
| 9 | tmax_avg_09 | ... September | 9 月上游最高月气温（平均） | average | ℃×10 |
| 10 | tmax_avg_10 | ... October | 10 月上游最高月气温（平均） | average | ℃×10 |
| 11 | tmax_avg_11 | ... November | 11 月上游最高月气温（平均） | average | ℃×10 |
| 12 | tmax_avg_12 | ... December | 12 月上游最高月气温（平均） | average | ℃×10 |

来源：WorldClim。

---

### 6) 月降水量（上游总和）— `monthly_prec_sum.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | prec_sum_01 | Sum of monthly precipitation for January | 1 月上游降水量总和 | sum | mm |
| 2 | prec_sum_02 | ... February | 2 月上游降水量总和 | sum | mm |
| 3 | prec_sum_03 | ... March | 3 月上游降水量总和 | sum | mm |
| 4 | prec_sum_04 | ... April | 4 月上游降水量总和 | sum | mm |
| 5 | prec_sum_05 | ... May | 5 月上游降水量总和 | sum | mm |
| 6 | prec_sum_06 | ... June | 6 月上游降水量总和 | sum | mm |
| 7 | prec_sum_07 | ... July | 7 月上游降水量总和 | sum | mm |
| 8 | prec_sum_08 | ... August | 8 月上游降水量总和 | sum | mm |
| 9 | prec_sum_09 | ... September | 9 月上游降水量总和 | sum | mm |
| 10 | prec_sum_10 | ... October | 10 月上游降水量总和 | sum | mm |
| 11 | prec_sum_11 | ... November | 11 月上游降水量总和 | sum | mm |
| 12 | prec_sum_12 | ... December | 12 月上游降水量总和 | sum | mm |

来源：WorldClim。

---

### 7) 月均最低气温（上游加权平均）— `monthly_tmin_weighted_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | tmin_wavg_01 | Min. monthly air temperature January | 1 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 2 | tmin_wavg_02 | ... February | 2 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 3 | tmin_wavg_03 | ... March | 3 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 4 | tmin_wavg_04 | ... April | 4 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 5 | tmin_wavg_05 | ... May | 5 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 6 | tmin_wavg_06 | ... June | 6 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 7 | tmin_wavg_07 | ... July | 7 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 8 | tmin_wavg_08 | ... August | 8 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 9 | tmin_wavg_09 | ... September | 9 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 10 | tmin_wavg_10 | ... October | 10 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 11 | tmin_wavg_11 | ... November | 11 月上游最低月气温（加权平均） | weighted average | ℃×10 |
| 12 | tmin_wavg_12 | ... December | 12 月上游最低月气温（加权平均） | weighted average | ℃×10 |

来源：WorldClim。

---

### 8) 月均最高气温（上游加权平均）— `monthly_tmax_weighted_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | tmax_wavg_01 | Max. monthly air temperature January | 1 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 2 | tmax_wavg_02 | ... February | 2 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 3 | tmax_wavg_03 | ... March | 3 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 4 | tmax_wavg_04 | ... April | 4 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 5 | tmax_wavg_05 | ... May | 5 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 6 | tmax_wavg_06 | ... June | 6 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 7 | tmax_wavg_07 | ... July | 7 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 8 | tmax_wavg_08 | ... August | 8 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 9 | tmax_wavg_09 | ... September | 9 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 10 | tmax_wavg_10 | ... October | 10 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 11 | tmax_wavg_11 | ... November | 11 月上游最高月气温（加权平均） | weighted average | ℃×10 |
| 12 | tmax_wavg_12 | ... December | 12 月上游最高月气温（加权平均） | weighted average | ℃×10 |

来源：WorldClim。

---

### 9) 月降水量（上游加权总和）— `monthly_prec_weighted_sum.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | prec_wsum_01 | Sum of monthly precipitation January | 1 月上游降水量加权总和 | weighted sum | mm |
| 2 | prec_wsum_02 | ... February | 2 月上游降水量加权总和 | weighted sum | mm |
| 3 | prec_wsum_03 | ... March | 3 月上游降水量加权总和 | weighted sum | mm |
| 4 | prec_wsum_04 | ... April | 4 月上游降水量加权总和 | weighted sum | mm |
| 5 | prec_wsum_05 | ... May | 5 月上游降水量加权总和 | weighted sum | mm |
| 6 | prec_wsum_06 | ... June | 6 月上游降水量加权总和 | weighted sum | mm |
| 7 | prec_wsum_07 | ... July | 7 月上游降水量加权总和 | weighted sum | mm |
| 8 | prec_wsum_08 | ... August | 8 月上游降水量加权总和 | weighted sum | mm |
| 9 | prec_wsum_09 | ... September | 9 月上游降水量加权总和 | weighted sum | mm |
| 10 | prec_wsum_10 | ... October | 10 月上游降水量加权总和 | weighted sum | mm |
| 11 | prec_wsum_11 | ... November | 11 月上游降水量加权总和 | weighted sum | mm |
| 12 | prec_wsum_12 | ... December | 12 月上游降水量加权总和 | weighted sum | mm |

来源：WorldClim。

---

### 10) 上游 Bioclim（平均/总和）— `hydroclim_average+sum.nc`

变量与 Bioclim 对应关系（见 ReadMe：hydro_*_01 ~ hydro_*_19）：

| Band | 变量代码 | 英文含义（Bioclim） | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | hydro_avg_01 | Annual Mean Upstream Temperature (BIO1) | 年均上游气温 | average | ℃×10 |
| 2 | hydro_avg_02 | Mean Upstream Diurnal Range (BIO2) | 上游昼夜温差均值（逐月最大-最小之均值） | average | ℃×10 |
| 3 | hydro_avg_03 | Upstream Isothermality (BIO3) | 上游等温性指数（02/07×100） | average | ×100 |
| 4 | hydro_avg_04 | Temperature Seasonality (BIO4) | 上游气温季节性（标准差×100） | average | ℃×10 |
| 5 | hydro_avg_05 | Max Temp of Warmest Month (BIO5) | 上游最暖月最高气温 | average | ℃×10 |
| 6 | hydro_avg_06 | Min Temp of Coldest Month (BIO6) | 上游最冷月最低气温 | average | ℃×10 |
| 7 | hydro_avg_07 | Temperature Annual Range (BIO7) | 上游年温差（BIO5-BIO6） | average | ℃×10 |
| 8 | hydro_avg_08 | Mean Temp of Wettest Quarter (BIO8) | 上游最湿季平均气温 | average | ℃×10 |
| 9 | hydro_avg_09 | Mean Temp of Driest Quarter (BIO9) | 上游最干季平均气温 | average | ℃×10 |
| 10 | hydro_avg_10 | Mean Temp of Warmest Quarter (BIO10) | 上游最暖季平均气温 | average | ℃×10 |
| 11 | hydro_avg_11 | Mean Temp of Coldest Quarter (BIO11) | 上游最冷季平均气温 | average | ℃×10 |
| 12 | hydro_avg_12 | Annual Upstream Precipitation (BIO12) | 年度上游降水量 | sum | mm |
| 13 | hydro_avg_13 | Precipitation of Wettest Month (BIO13) | 上游最湿月降水量 | sum | mm |
| 14 | hydro_avg_14 | Precipitation of Driest Month (BIO14) | 上游最干月降水量 | sum | mm |
| 15 | hydro_avg_15 | Precipitation Seasonality (BIO15) | 上游降水季节性（变异系数） | sum | ×100 |
| 16 | hydro_avg_16 | Precipitation of Wettest Quarter (BIO16) | 上游最湿季降水量 | sum | mm |
| 17 | hydro_avg_17 | Precipitation of Driest Quarter (BIO17) | 上游最干季降水量 | sum | mm |
| 18 | hydro_avg_18 | Precipitation of Warmest Quarter (BIO18) | 上游最暖季降水量 | sum | mm |
| 19 | hydro_avg_19 | Precipitation of Coldest Quarter (BIO19) | 上游最冷季降水量 | sum | mm |

来源：WorldClim。

---

### 11) 上游 Bioclim（加权平均/总和）— `hydroclim_weighted_average+sum.nc`

| Band | 变量代码 | 英文含义（Bioclim） | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | hydro_wavg_01 | Annual Mean Upstream Temperature (BIO1) | 年均上游气温（加权） | weighted average | ℃×10 |
| 2 | hydro_wavg_02 | Mean Upstream Diurnal Range (BIO2) | 上游昼夜温差均值（加权） | weighted average | ℃×10 |
| 3 | hydro_wavg_03 | Upstream Isothermality (BIO3) | 上游等温性指数（加权；02/07×100） | weighted average | ×100 |
| 4 | hydro_wavg_04 | Temperature Seasonality (BIO4) | 上游气温季节性（加权；标准差×100） | weighted average | ℃×10 |
| 5 | hydro_wavg_05 | Max Temp of Warmest Month (BIO5) | 上游最暖月最高气温（加权） | weighted average | ℃×10 |
| 6 | hydro_wavg_06 | Min Temp of Coldest Month (BIO6) | 上游最冷月最低气温（加权） | weighted average | ℃×10 |
| 7 | hydro_wavg_07 | Temperature Annual Range (BIO7) | 上游年温差（加权；BIO5-BIO6） | weighted average | ℃×10 |
| 8 | hydro_wavg_08 | Mean Temp of Wettest Quarter (BIO8) | 上游最湿季平均气温（加权） | weighted average | ℃×10 |
| 9 | hydro_wavg_09 | Mean Temp of Driest Quarter (BIO9) | 上游最干季平均气温（加权） | weighted average | ℃×10 |
| 10 | hydro_wavg_10 | Mean Temp of Warmest Quarter (BIO10) | 上游最暖季平均气温（加权） | weighted average | ℃×10 |
| 11 | hydro_wavg_11 | Mean Temp of Coldest Quarter (BIO11) | 上游最冷季平均气温（加权） | weighted average | ℃×10 |
| 12 | hydro_wavg_12 | Annual Upstream Precipitation (BIO12) | 年度上游降水量（加权） | weighted sum | mm |
| 13 | hydro_wavg_13 | Precipitation of Wettest Month (BIO13) | 上游最湿月降水量（加权） | weighted sum | mm |
| 14 | hydro_wavg_14 | Precipitation of Driest Month (BIO14) | 上游最干月降水量（加权） | weighted sum | mm |
| 15 | hydro_wavg_15 | Precipitation Seasonality (BIO15) | 上游降水季节性（加权；变异系数） | weighted sum | ×100 |
| 16 | hydro_wavg_16 | Precipitation of Wettest Quarter (BIO16) | 上游最湿季降水量（加权） | weighted sum | mm |
| 17 | hydro_wavg_17 | Precipitation of Driest Quarter (BIO17) | 上游最干季降水量（加权） | weighted sum | mm |
| 18 | hydro_wavg_18 | Precipitation of Warmest Quarter (BIO18) | 上游最暖季降水量（加权） | weighted sum | mm |
| 19 | hydro_wavg_19 | Precipitation of Coldest Quarter (BIO19) | 上游最冷季降水量（加权） | weighted sum | mm |

来源：WorldClim。

---

### 12) 土地覆盖（最小）— `landcover_minimum.nc`

注：各类为 Consensus Land Cover（CLC）12 类。最小值表示上游区域内该类所占比例的最小统计。

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | lc_min_01 | Evergreen/deciduous needleleaf trees | 常绿/落叶针叶林（最小） | minimum | % |
| 2 | lc_min_02 | Evergreen broadleaf trees | 常绿阔叶林（最小） | minimum | % |
| 3 | lc_min_03 | Deciduous broadleaf trees | 落叶阔叶林（最小） | minimum | % |
| 4 | lc_min_04 | Mixed/other trees | 混交/其他树林（最小） | minimum | % |
| 5 | lc_min_05 | Shrubs | 灌木（最小） | minimum | % |
| 6 | lc_min_06 | Herbaceous vegetation | 草本植被（最小） | minimum | % |
| 7 | lc_min_07 | Cultivated and managed vegetation | 农业耕作与管理植被（最小） | minimum | % |
| 8 | lc_min_08 | Regularly flooded shrub/herbaceous | 经常性淹涝灌/草植被（最小） | minimum | % |
| 9 | lc_min_09 | Urban/built-up | 城市/建成区（最小） | minimum | % |
| 10 | lc_min_10 | Snow/ice | 积雪/冰（最小） | minimum | % |
| 11 | lc_min_11 | Barren lands/sparse vegetation | 荒漠/稀疏植被（最小） | minimum | % |
| 12 | lc_min_12 | Open water | 开阔水体（最小） | minimum | % |

来源：CLC。

---

### 13) 土地覆盖（最大）— `landcover_maximum.nc`

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | lc_max_01 | Evergreen/deciduous needleleaf trees | 常绿/落叶针叶林（最大） | maximum | % |
| 2 | lc_max_02 | Evergreen broadleaf trees | 常绿阔叶林（最大） | maximum | % |
| 3 | lc_max_03 | Deciduous broadleaf trees | 落叶阔叶林（最大） | maximum | % |
| 4 | lc_max_04 | Mixed/other trees | 混交/其他树林（最大） | maximum | % |
| 5 | lc_max_05 | Shrubs | 灌木（最大） | maximum | % |
| 6 | lc_max_06 | Herbaceous vegetation | 草本植被（最大） | maximum | % |
| 7 | lc_max_07 | Cultivated and managed vegetation | 农业耕作与管理植被（最大） | maximum | % |
| 8 | lc_max_08 | Regularly flooded shrub/herbaceous | 经常性淹涝灌/草植被（最大） | maximum | % |
| 9 | lc_max_09 | Urban/built-up | 城市/建成区（最大） | maximum | % |
| 10 | lc_max_10 | Snow/ice | 积雪/冰（最大） | maximum | % |
| 11 | lc_max_11 | Barren lands/sparse vegetation | 荒漠/稀疏植被（最大） | maximum | % |
| 12 | lc_max_12 | Open water | 开阔水体（最大） | maximum | % |

来源：CLC。

---

### 14) 土地覆盖（范围）— `landcover_range.nc`

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | lc_range_01 | Evergreen/deciduous needleleaf trees | 常绿/落叶针叶林（范围） | range | % |
| 2 | lc_range_02 | Evergreen broadleaf trees | 常绿阔叶林（范围） | range | % |
| 3 | lc_range_03 | Deciduous broadleaf trees | 落叶阔叶林（范围） | range | % |
| 4 | lc_range_04 | Mixed/other trees | 混交/其他树林（范围） | range | % |
| 5 | lc_range_05 | Shrubs | 灌木（范围） | range | % |
| 6 | lc_range_06 | Herbaceous vegetation | 草本植被（范围） | range | % |
| 7 | lc_range_07 | Cultivated and managed vegetation | 农业耕作与管理植被（范围） | range | % |
| 8 | lc_range_08 | Regularly flooded shrub/herbaceous | 经常性淹涝灌/草植被（范围） | range | % |
| 9 | lc_range_09 | Urban/built-up | 城市/建成区（范围） | range | % |
| 10 | lc_range_10 | Snow/ice | 积雪/冰（范围） | range | % |
| 11 | lc_range_11 | Barren lands/sparse vegetation | 荒漠/稀疏植被（范围） | range | % |
| 12 | lc_range_12 | Open water | 开阔水体（范围） | range | % |

来源：CLC。

---

### 15) 土地覆盖（平均）— `landcover_average.nc`

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | lc_avg_01 | Evergreen/deciduous needleleaf trees | 常绿/落叶针叶林（平均） | average | % |
| 2 | lc_avg_02 | Evergreen broadleaf trees | 常绿阔叶林（平均） | average | % |
| 3 | lc_avg_03 | Deciduous broadleaf trees | 落叶阔叶林（平均） | average | % |
| 4 | lc_avg_04 | Mixed/other trees | 混交/其他树林（平均） | average | % |
| 5 | lc_avg_05 | Shrubs | 灌木（平均） | average | % |
| 6 | lc_avg_06 | Herbaceous vegetation | 草本植被（平均） | average | % |
| 7 | lc_avg_07 | Cultivated and managed vegetation | 农业耕作与管理植被（平均） | average | % |
| 8 | lc_avg_08 | Regularly flooded shrub/herbaceous | 经常性淹涝灌/草植被（平均） | average | % |
| 9 | lc_avg_09 | Urban/built-up | 城市/建成区（平均） | average | % |
| 10 | lc_avg_10 | Snow/ice | 积雪/冰（平均） | average | % |
| 11 | lc_avg_11 | Barren lands/sparse vegetation | 荒漠/稀疏植被（平均） | average | % |
| 12 | lc_avg_12 | Open water | 开阔水体（平均） | average | % |

来源：CLC。

---

### 16) 土地覆盖（加权平均）— `landcover_weighted_average.nc`

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | lc_wavg_01 | Evergreen/deciduous needleleaf trees | 常绿/落叶针叶林（加权平均） | weighted average | % |
| 2 | lc_wavg_02 | Evergreen broadleaf trees | 常绿阔叶林（加权平均） | weighted average | % |
| 3 | lc_wavg_03 | Deciduous broadleaf trees | 落叶阔叶林（加权平均） | weighted average | % |
| 4 | lc_wavg_04 | Mixed/other trees | 混交/其他树林（加权平均） | weighted average | % |
| 5 | lc_wavg_05 | Shrubs | 灌木（加权平均） | weighted average | % |
| 6 | lc_wavg_06 | Herbaceous vegetation | 草本植被（加权平均） | weighted average | % |
| 7 | lc_wavg_07 | Cultivated and managed vegetation | 农业耕作与管理植被（加权平均） | weighted average | % |
| 8 | lc_wavg_08 | Regularly flooded shrub/herbaceous | 经常性淹涝灌/草植被（加权平均） | weighted average | % |
| 9 | lc_wavg_09 | Urban/built-up | 城市/建成区（加权平均） | weighted average | % |
| 10 | lc_wavg_10 | Snow/ice | 积雪/冰（加权平均） | weighted average | % |
| 11 | lc_wavg_11 | Barren lands/sparse vegetation | 荒漠/稀疏植被（加权平均） | weighted average | % |
| 12 | lc_wavg_12 | Open water | 开阔水体（加权平均） | weighted average | % |

来源：CLC。

---

### 17) 地表地质（上游加权计数）— `geology_weighted_sum.nc`

注：表示上游范围内各地质年代/类型的加权栅格计数（weighted count）。中文对照采用标准地质年代名称。

| Band | 变量代码 | 英文类名 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | geo_wsum_01 | Archean | 太古代 | weighted sum | 加权计数 |
| 2 | geo_wsum_02 | Archean, Permian | 太古代 + 二叠纪 | weighted sum | 加权计数 |
| 3 | geo_wsum_03 | Cambrian | 寒武纪 | weighted sum | 加权计数 |
| 4 | geo_wsum_04 | Cambrian, Precambrian | 寒武纪 + 前寒武纪 | weighted sum | 加权计数 |
| 5 | geo_wsum_05 | Cambrian, Proterozoic | 寒武纪 + 元古代 | weighted sum | 加权计数 |
| 6 | geo_wsum_06 | Carboniferous | 石炭纪 | weighted sum | 加权计数 |
| 7 | geo_wsum_07 | Carboniferous, Devonian | 石炭纪 + 泥盆纪 | weighted sum | 加权计数 |
| 8 | geo_wsum_08 | Carboniferous, Miocene | 石炭纪 + 中新世 | weighted sum | 加权计数 |
| 9 | geo_wsum_09 | Cenozoic | 新生代 | weighted sum | 加权计数 |
| 10 | geo_wsum_10 | Cenozoic, Mesozoic | 新生代 + 中生代 | weighted sum | 加权计数 |
| 11 | geo_wsum_11 | Cretaceous | 白垩纪 | weighted sum | 加权计数 |
| 12 | geo_wsum_12 | Cretaceous, Carboniferous | 白垩纪 + 石炭纪 | weighted sum | 加权计数 |
| 13 | geo_wsum_13 | Cretaceous, Devonian | 白垩纪 + 泥盆纪 | weighted sum | 加权计数 |
| 14 | geo_wsum_14 | Cretaceous, Jurassic | 白垩纪 + 侏罗纪 | weighted sum | 加权计数 |
| 15 | geo_wsum_15 | Cretaceous, Mississippian | 白垩纪 + 密西西比世 | weighted sum | 加权计数 |
| 16 | geo_wsum_16 | Cretaceous, Paleogene, Neogene | 白垩纪 + 古近纪 + 新近纪 | weighted sum | 加权计数 |
| 17 | geo_wsum_17 | Cretaceous, Permian | 白垩纪 + 二叠纪 | weighted sum | 加权计数 |
| 18 | geo_wsum_18 | Cretaceous, Tertiary | 白垩纪 + 第三纪 | weighted sum | 加权计数 |
| 19 | geo_wsum_19 | Cretaceous, Triassic | 白垩纪 + 三叠纪 | weighted sum | 加权计数 |
| 20 | geo_wsum_20 | Devonian | 泥盆纪 | weighted sum | 加权计数 |
| 21 | geo_wsum_21 | Devonian, Cambrian | 泥盆纪 + 寒武纪 | weighted sum | 加权计数 |
| 22 | geo_wsum_22 | Devonian, Ordovician | 泥盆纪 + 奥陶纪 | weighted sum | 加权计数 |
| 23 | geo_wsum_23 | Devonian, Proterozoic | 泥盆纪 + 元古代 | weighted sum | 加权计数 |
| 24 | geo_wsum_24 | Devonian, Silurian | 泥盆纪 + 志留纪 | weighted sum | 加权计数 |
| 25 | geo_wsum_25 | Devonian, Silurian, Ordovician | 泥盆纪 + 志留纪 + 奥陶纪 | weighted sum | 加权计数 |
| 26 | geo_wsum_26 | Holocene | 全新世 | weighted sum | 加权计数 |
| 27 | geo_wsum_27 | Ice | 冰 | weighted sum | 加权计数 |
| 28 | geo_wsum_28 | Jurassic | 侏罗纪 | weighted sum | 加权计数 |
| 29 | geo_wsum_29 | Jurassic, Cambrian | 侏罗纪 + 寒武纪 | weighted sum | 加权计数 |
| 30 | geo_wsum_30 | Jurassic, Carboniferous | 侏罗纪 + 石炭纪 | weighted sum | 加权计数 |
| 31 | geo_wsum_31 | Jurassic, Devonian | 侏罗纪 + 泥盆纪 | weighted sum | 加权计数 |
| 32 | geo_wsum_32 | Jurassic, Mississippian | 侏罗纪 + 密西西比世 | weighted sum | 加权计数 |
| 33 | geo_wsum_33 | Jurassic, Ordovician | 侏罗纪 + 奥陶纪 | weighted sum | 加权计数 |
| 34 | geo_wsum_34 | Jurassic, Permian | 侏罗纪 + 二叠纪 | weighted sum | 加权计数 |
| 35 | geo_wsum_35 | Jurassic, Triassic | 侏罗纪 + 三叠纪 | weighted sum | 加权计数 |
| 36 | geo_wsum_36 | Kimberlite | 金伯利岩 | weighted sum | 加权计数 |
| 37 | geo_wsum_37 | Mesozoic | 中生代 | weighted sum | 加权计数 |
| 38 | geo_wsum_38 | Mesozoic, Cenozoic | 中生代 + 新生代 | weighted sum | 加权计数 |
| 39 | geo_wsum_39 | Mesozoic, Paleozoic | 中生代 + 古生代 | weighted sum | 加权计数 |
| 40 | geo_wsum_40 | Mesozoic, Paleozoic | 中生代 + 古生代 | weighted sum | 加权计数 |
| 41 | geo_wsum_41 | Miocene | 中新世 | weighted sum | 加权计数 |
| 42 | geo_wsum_42 | Mississippian | 密西西比世（早石炭世） | weighted sum | 加权计数 |
| 43 | geo_wsum_43 | Mississippian, Cambrian | 密西西比世 + 寒武纪 | weighted sum | 加权计数 |
| 44 | geo_wsum_44 | Mississippian, Devonian | 密西西比世 + 泥盆纪 | weighted sum | 加权计数 |
| 45 | geo_wsum_45 | Neogene | 新近纪 | weighted sum | 加权计数 |
| 46 | geo_wsum_46 | Neogene, Paleogene | 新近纪 + 古近纪 | weighted sum | 加权计数 |
| 47 | geo_wsum_47 | Ordovician | 奥陶纪 | weighted sum | 加权计数 |
| 48 | geo_wsum_48 | Ordovician, Cambrian | 奥陶纪 + 寒武纪 | weighted sum | 加权计数 |
| 49 | geo_wsum_49 | Ordovician, Proterozoic | 奥陶纪 + 元古代 | weighted sum | 加权计数 |
| 50 | geo_wsum_50 | Paleogene | 古近纪 | weighted sum | 加权计数 |
| 51 | geo_wsum_51 | Paleogene, Cretaceous | 古近纪 + 白垩纪 | weighted sum | 加权计数 |
| 52 | geo_wsum_52 | Paleozoic | 古生代 | weighted sum | 加权计数 |
| 53 | geo_wsum_53 | Paleozoic, Mesozoic | 古生代 + 中生代 | weighted sum | 加权计数 |
| 54 | geo_wsum_54 | Paleozoic, Precambrian | 古生代 + 前寒武纪 | weighted sum | 加权计数 |
| 55 | geo_wsum_55 | Paleozoic, Proterozoic | 古生代 + 元古代 | weighted sum | 加权计数 |
| 56 | geo_wsum_56 | Pennsylvanian | 宾夕法尼亚世（晚石炭世） | weighted sum | 加权计数 |
| 57 | geo_wsum_57 | Pennsylvanian, Devonian | 宾夕法尼亚世 + 泥盆纪 | weighted sum | 加权计数 |
| 58 | geo_wsum_58 | Pennsylvanian, Mississippian | 宾夕法尼亚世 + 密西西比世 | weighted sum | 加权计数 |
| 59 | geo_wsum_59 | Permian | 二叠纪 | weighted sum | 加权计数 |
| 60 | geo_wsum_60 | Permian, Carboniferous | 二叠纪 + 石炭纪 | weighted sum | 加权计数 |
| 61 | geo_wsum_61 | Permian, Devonian | 二叠纪 + 泥盆纪 | weighted sum | 加权计数 |
| 62 | geo_wsum_62 | Permian, Mississippian | 二叠纪 + 密西西比世 | weighted sum | 加权计数 |
| 63 | geo_wsum_63 | Permian, Pennsylvanian | 二叠纪 + 宾夕法尼亚世 | weighted sum | 加权计数 |
| 64 | geo_wsum_64 | Permian, Triassic | 二叠纪 + 三叠纪 | weighted sum | 加权计数 |
| 65 | geo_wsum_65 | Pleistocene | 更新世 | weighted sum | 加权计数 |
| 66 | geo_wsum_66 | Pliocene, Quaternary | 上新世 + 第四纪 | weighted sum | 加权计数 |
| 67 | geo_wsum_67 | Precambrian | 前寒武纪 | weighted sum | 加权计数 |
| 68 | geo_wsum_68 | Precambrian, Devonian | 前寒武纪 + 泥盆纪 | weighted sum | 加权计数 |
| 69 | geo_wsum_69 | Precambrian, Paleozoic | 前寒武纪 + 古生代 | weighted sum | 加权计数 |
| 70 | geo_wsum_70 | Proterozoic | 元古代 | weighted sum | 加权计数 |
| 71 | geo_wsum_71 | Proterozoic, Archean | 元古代 + 太古代 | weighted sum | 加权计数 |
| 72 | geo_wsum_72 | Quaternary | 第四纪 | weighted sum | 加权计数 |
| 73 | geo_wsum_73 | Quaternary, Neogene | 第四纪 + 新近纪 | weighted sum | 加权计数 |
| 74 | geo_wsum_74 | Quaternary, Tertiary | 第四纪 + 第三纪 | weighted sum | 加权计数 |
| 75 | geo_wsum_75 | Salt | 盐类 | weighted sum | 加权计数 |
| 76 | geo_wsum_76 | Silurian | 志留纪 | weighted sum | 加权计数 |
| 77 | geo_wsum_77 | Silurian, Cambrian | 志留纪 + 寒武纪 | weighted sum | 加权计数 |
| 78 | geo_wsum_78 | Silurian, Ordovician | 志留纪 + 奥陶纪 | weighted sum | 加权计数 |
| 79 | geo_wsum_79 | Silurian, Proterozoic | 志留纪 + 元古代 | weighted sum | 加权计数 |
| 80 | geo_wsum_80 | Tertiary | 第三纪 | weighted sum | 加权计数 |
| 81 | geo_wsum_81 | Tertiary, Cretaceous | 第三纪 + 白垩纪 | weighted sum | 加权计数 |
| 82 | geo_wsum_82 | Triassic | 三叠纪 | weighted sum | 加权计数 |
| 83 | geo_wsum_83 | Triassic, Carboniferous | 三叠纪 + 石炭纪 | weighted sum | 加权计数 |
| 84 | geo_wsum_84 | Triassic, Devonian | 三叠纪 + 泥盆纪 | weighted sum | 加权计数 |
| 85 | geo_wsum_85 | Triassic, Mississippian | 三叠纪 + 密西西比世 | weighted sum | 加权计数 |
| 86 | geo_wsum_86 | Triassic, Ordovician | 三叠纪 + 奥陶纪 | weighted sum | 加权计数 |
| 87 | geo_wsum_87 | Triassic, Paleozoic | 三叠纪 + 古生代 | weighted sum | 加权计数 |
| 88 | geo_wsum_88 | Triassic, Pennsylvanian | 三叠纪 + 宾夕法尼亚世 | weighted sum | 加权计数 |
| 89 | geo_wsum_89 | Triassic, Permian | 三叠纪 + 二叠纪 | weighted sum | 加权计数 |
| 90 | geo_wsum_90 | Triassic, Proterozoic | 三叠纪 + 元古代 | weighted sum | 加权计数 |
| 91 | geo_wsum_91 | Unknown | 未知 | weighted sum | 加权计数 |
| 92 | geo_wsum_92 | Water | 水体 | weighted sum | 加权计数 |

来源：USGS 地表地质。

---

### 18) 土壤（最小）— `soil_minimum.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | soil_min_01 | Soil organic carbon | 土壤有机碳（最小） | minimum | g/kg |
| 2 | soil_min_02 | Soil pH in H2O | 土壤 pH（最小，水浸提；×10） | minimum | pH×10 |
| 3 | soil_min_03 | Sand content mass fraction | 土壤砂粒含量（最小） | minimum | % |
| 4 | soil_min_04 | Silt content mass fraction | 土壤粉粒含量（最小） | minimum | % |
| 5 | soil_min_05 | Clay content mass fraction | 土壤黏粒含量（最小） | minimum | % |
| 6 | soil_min_06 | Coarse fragments (>2 mm) volumetric | 大颗粒体积分数>2mm（最小） | minimum | % |
| 7 | soil_min_07 | Cation exchange capacity | 阳离子交换量（最小） | minimum | cmol/kg |
| 8 | soil_min_08 | Bulk density of fine earth | 细土体积密度（最小） | minimum | kg/m³ |
| 9 | soil_min_09 | Depth to bedrock (≤240 cm) | 基岩深度（最小，≤240 cm） | minimum | cm |
| 10 | soil_min_10 | Probability of R horizon | 基岩（R 层）出现概率（最小） | minimum | % |

来源：ISRIC SoilGrids。

---

### 19) 土壤（最大）— `soil_maximum.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | soil_max_01 | Soil organic carbon | 土壤有机碳（最大） | maximum | g/kg |
| 2 | soil_max_02 | Soil pH in H2O | 土壤 pH（最大，水浸提；×10） | maximum | pH×10 |
| 3 | soil_max_03 | Sand content mass fraction | 土壤砂粒含量（最大） | maximum | % |
| 4 | soil_max_04 | Silt content mass fraction | 土壤粉粒含量（最大） | maximum | % |
| 5 | soil_max_05 | Clay content mass fraction | 土壤黏粒含量（最大） | maximum | % |
| 6 | soil_max_06 | Coarse fragments (>2 mm) volumetric | 大颗粒体积分数>2mm（最大） | maximum | % |
| 7 | soil_max_07 | Cation exchange capacity | 阳离子交换量（最大） | maximum | cmol/kg |
| 8 | soil_max_08 | Bulk density of fine earth | 细土体积密度（最大） | maximum | kg/m³ |
| 9 | soil_max_09 | Depth to bedrock (≤240 cm) | 基岩深度（最大，≤240 cm） | maximum | cm |
| 10 | soil_max_10 | Probability of R horizon | 基岩（R 层）出现概率（最大） | maximum | % |

来源：ISRIC SoilGrids。

---

### 20) 土壤（范围）— `soil_range.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | soil_range_01 | Soil organic carbon | 土壤有机碳（范围） | range | g/kg |
| 2 | soil_range_02 | Soil pH in H2O | 土壤 pH（范围，×10） | range | pH×10 |
| 3 | soil_range_03 | Sand content mass fraction | 土壤砂粒含量（范围） | range | % |
| 4 | soil_range_04 | Silt content mass fraction | 土壤粉粒含量（范围） | range | % |
| 5 | soil_range_05 | Clay content mass fraction | 土壤黏粒含量（范围） | range | % |
| 6 | soil_range_06 | Coarse fragments (>2 mm) volumetric | 大颗粒体积分数>2mm（范围） | range | % |
| 7 | soil_range_07 | Cation exchange capacity | 阳离子交换量（范围） | range | cmol/kg |
| 8 | soil_range_08 | Bulk density of fine earth | 细土体积密度（范围） | range | kg/m³ |
| 9 | soil_range_09 | Depth to bedrock (≤240 cm) | 基岩深度（范围，≤240 cm） | range | cm |
| 10 | soil_range_10 | Probability of R horizon | 基岩（R 层）出现概率（范围） | range | % |

来源：ISRIC SoilGrids。

---

### 21) 土壤（平均）— `soil_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | soil_avg_01 | Soil organic carbon | 土壤有机碳（平均） | average | g/kg |
| 2 | soil_avg_02 | Soil pH in H2O | 土壤 pH（平均，×10） | average | pH×10 |
| 3 | soil_avg_03 | Sand content mass fraction | 土壤砂粒含量（平均） | average | % |
| 4 | soil_avg_04 | Silt content mass fraction | 土壤粉粒含量（平均） | average | % |
| 5 | soil_avg_05 | Clay content mass fraction | 土壤黏粒含量（平均） | average | % |
| 6 | soil_avg_06 | Coarse fragments (>2 mm) volumetric | 大颗粒体积分数>2mm（平均） | average | % |
| 7 | soil_avg_07 | Cation exchange capacity | 阳离子交换量（平均） | average | cmol/kg |
| 8 | soil_avg_08 | Bulk density of fine earth | 细土体积密度（平均） | average | kg/m³ |
| 9 | soil_avg_09 | Depth to bedrock (≤240 cm) | 基岩深度（平均，≤240 cm） | average | cm |
| 10 | soil_avg_10 | Probability of R horizon | 基岩（R 层）出现概率（平均） | average | % |

来源：ISRIC SoilGrids。

---

### 22) 土壤（加权平均）— `soil_weighted_average.nc`

| Band | 变量代码 | 英文含义 | 中文解释 | 聚合 | 单位 |
|---:|---|---|---|---|---|
| 1 | soil_wavg_01 | Soil organic carbon | 土壤有机碳（加权平均） | weighted average | g/kg |
| 2 | soil_wavg_02 | Soil pH in H2O | 土壤 pH（加权平均，×10） | weighted average | pH×10 |
| 3 | soil_wavg_03 | Sand content mass fraction | 土壤砂粒含量（加权平均） | weighted average | % |
| 4 | soil_wavg_04 | Silt content mass fraction | 土壤粉粒含量（加权平均） | weighted average | % |
| 5 | soil_wavg_05 | Clay content mass fraction | 土壤黏粒含量（加权平均） | weighted average | % |
| 6 | soil_wavg_06 | Coarse fragments (>2 mm) volumetric | 大颗粒体积分数>2mm（加权平均） | weighted average | % |
| 7 | soil_wavg_07 | Cation exchange capacity | 阳离子交换量（加权平均） | weighted average | cmol/kg |
| 8 | soil_wavg_08 | Bulk density of fine earth | 细土体积密度（加权平均） | weighted average | kg/m³ |
| 9 | soil_wavg_09 | Depth to bedrock (≤240 cm) | 基岩深度（加权平均，≤240 cm） | weighted average | cm |
| 10 | soil_wavg_10 | Probability of R horizon | 基岩（R 层）出现概率（加权平均） | weighted average | % |

来源：ISRIC SoilGrids。

---

### 23) 质量控制 — `quality_control.nc`

| Band | 变量代码 | 英文含义 | 中文解释 |
|---:|---|---|---|
| 1 | missing_cells | Cells filled with max. neighbour value | 采用邻域最大值填补的像元 |
| 2 | cells_removed | Cells that were removed manually | 被人工剔除的像元 |

---

以上共计 324 个变量。若需批量筛选变量或导出变量清单为 CSV，请告知需要的筛选规则（例如“仅温度相关”“仅加权版本”“去除土地覆盖”），我可为您自动生成对应的提取与导出脚本（R 或 Python），并将结果保存到 `output` 目录以便后续制图与分析（保证 Nature 期刊标准的数据与图件规范）。


