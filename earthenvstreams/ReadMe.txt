
Near-global freshwater-specific environmental variables for biodiversity analyses in 1km resolution 
====================================================================================================

version 1.0, October 2015

Sami Domisch, Giuseppe Amatulli, Walter Jetz
Department of Ecology and Evolutionary Biology, Yale University, 165 Prospect Street, New Haven, CT 06511, USA
Corresponding author: Sami Domisch (sami.domisch@yale.edu)



General information:
======================
All 324 layers are available as near-global netCDF-4 files, where each file contains one 
variable and aggregation type (e.g. all 12 average landcover layers). See the example R-script 
in the Supplementary Information to load the variables into R.

Extent: (-145, -56, 180, 60)
Cell size: 0.0083333333 degree (30 arcsec)
Projection: WGS84
Number of columns: 39000
Number of rows: 13920
Pixel type: Byte, Int32, Float64
NoData value: -127 (Byte), -999 (Int32) or -9999 (Float64)


Citation:
=========
Domisch, S., Amatulli, G., Jetz, W. (2015). Near-global freshwater-specific environmental 
		variables for biodiversity analyses in 1km resolution. Scientific Data

		
Please also cite the HydroSHEDS and the GLWD data sets, and depending on the variables that are used, 
also the primary data source of these:

HydroSHEDS - Lehner, B., Verdin, K., & Jarvis, A. (2008). New global hydrography derived from 
		spaceborne elevation data. EOS, Transactions American Geophysical Union, 89(10), 93-94.

GLWD - Lehner, B., & Döll, P. (2004). Development and validation of a global database of lakes, 
		reservoirs and wetlands. Journal of Hydrology, 296(1), 1-22.

WorldClim - Hijmans, R. J., Cameron, S. E., Parra, J. L., Jones, P. G., & Jarvis, A. (2005). Very 
		high resolution interpolated climate surfaces for global land areas. International Journal 
		of Climatology, 25(15), 1965-1978.

Consensus Land cover (CLC) - Tuanmu, M. N., & Jetz, W. (2014). A global 1-km consensus land-cover product 
		for biodiversity and ecosystem modelling. Global Ecology and Biogeography, 23(9), 1031-1045.

Surface geology - USGS International Surface Geology < http://certmapper.cr.usgs.gov/data/envision/index.html?widgets=geologymaps >. 
		Downloaded from < http://worldgrids.org/doku.php?id=wiki:geaisg3 >

Soil - Hengl, T., de Jesus, J. M., MacMillan, R. A., Batjes, N.H., Heuvelink, G.B., Ribeiro, E., 
		Samuel-Rosa, A., Kempen, B., Leenaars, J. G. B., Walsh, M. G., Gonzalez & Gonzalez, M. R. (2014). 
		SoilGrids1km—Global soil information based on automated mapping. PLoS ONE, e105992




Hydroclimatic variables:
========================

The hydroclimatic variables ("hydro_*.nc") were created following the "bioclim" framework (see http://worldclim.org/bioclim for details), 
using the upstream temperature (average or weighted average) and upstream precipitation (sum or weighted sum):

hydro_*_01 = Annual Mean Upstream Temperature
hydro_*_02 = Mean Upstream Diurnal Range (Mean of monthly (max temp - min temp))
hydro_*_03 = Upstream Isothermality (hydro_02 / hydro_07) (* 100)
hydro_*_04 = Upstream Temperature Seasonality (standard deviation *100)
hydro_*_05 = Maximum Upstream Temperature of Warmest Month
hydro_*_06 = Minimum Upstream Temperature of Coldest Month
hydro_*_07 = Upstream Temperature Annual Range (hydro_05 - hydro_06)
hydro_*_08 = Mean Upstream Temperature of Wettest Quarter
hydro_*_09 = Mean Upstream Temperature of Driest Quarter
hydro_*_10 = Mean Upstream Temperature of Warmest Quarter
hydro_*_11 = Mean Upstream Temperature of Coldest Quarter
hydro_*_12 = Annual Upstream Precipitation
hydro_*_13 = Upstream Precipitation of Wettest Month
hydro_*_14 = Upstream Precipitation of Driest Month
hydro_*_15 = Upstream Precipitation Seasonality (Coefficient of Variation)
hydro_*_16 = Upstream Precipitation of Wettest Quarter
hydro_*_17 = Upstream Precipitation of Driest Quarter
hydro_*_18 = Upstream Precipitation of Warmest Quarter
hydro_*_19 = Upstream Precipitation of Coldest Quarter





Layer description:
==================



|netCDF file / Variable name       |Band ID |Variable code |Variable explanation                              |Aggregation      |Unit           |Type    |NoData |Source     |
|==================================|========|==============|==================================================|=================|===============|========|=======|===========|
|elevation.nc                      |1       |dem_min       |Minimum elevation                                 |minimum          |[m]            |Int32   |-999   |HydroSHEDS |
|elevation.nc                      |2       |dem_max       |Maximum elevation                                 |maximum          |[m]            |Int32   |-999   |HydroSHEDS |
|elevation.nc                      |3       |dem_range     |Elevation range                                   |range            |[m]            |Int32   |-999   |HydroSHEDS |
|elevation.nc                      |4       |dem_avg       |Average elevation                                 |average          |[m]            |Int32   |-999   |HydroSHEDS |
|slope.nc                          |1       |slope_min     |Minimum slope                                     |minimum          |[°] * 100      |Int32   |-999   |HydroSHEDS |
|slope.nc                          |2       |slope_max     |Maximum slope                                     |maximum          |[°] * 100      |Int32   |-999   |HydroSHEDS |
|slope.nc                          |3       |slope_range   |Slope range                                       |range            |[°] * 100      |Int32   |-999   |HydroSHEDS |
|slope.nc                          |4       |slope_avg     |Average slope                                     |average          |[°] * 100      |Int32   |-999   |HydroSHEDS |
|flow_acc.nc                       |1       |flow_length   |Number of upstream stream grid cells              |sum              |count          |Int32   |-999   |HydroSHEDS |
|flow_acc.nc                       |2       |flow_acc      |Number of upstream catchment grid cells           |sum              |count          |Int32   |-999   |HydroSHEDS |
|monthly_tmin_average.nc           |1       |tmin_avg_01   |Min. monthly air temperature for January          |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |2       |tmin_avg_02   |Min. monthly air temperature for February         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |3       |tmin_avg_03   |Min. monthly air temperature for March            |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |4       |tmin_avg_04   |Min. monthly air temperature for April            |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |5       |tmin_avg_05   |Min. monthly air temperature for May              |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |6       |tmin_avg_06   |Min. monthly air temperature for June             |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |7       |tmin_avg_07   |Min. monthly air temperature for July             |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |8       |tmin_avg_08   |Min. monthly air temperature for August           |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |9       |tmin_avg_09   |Min. monthly air temperature for September        |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |10      |tmin_avg_10   |Min. monthly air temperature for October          |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |11      |tmin_avg_11   |Min. monthly air temperature for November         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_average.nc           |12      |tmin_avg_12   |Min. monthly air temperature for December         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |1       |tmax_avg_01   |Max. monthly air temperature for January          |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |2       |tmax_avg_02   |Max. monthly air temperature for February         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |3       |tmax_avg_03   |Max. monthly air temperature for March            |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |4       |tmax_avg_04   |Max. monthly air temperature for April            |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |5       |tmax_avg_05   |Max. monthly air temperature for May              |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |6       |tmax_avg_06   |Max. monthly air temperature for June             |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |7       |tmax_avg_07   |Max. monthly air temperature for July             |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |8       |tmax_avg_08   |Max. monthly air temperature for August           |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |9       |tmax_avg_09   |Max. monthly air temperature for September        |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |10      |tmax_avg_10   |Max. monthly air temperature for October          |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |11      |tmax_avg_11   |Max. monthly air temperature for November         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_average.nc           |12      |tmax_avg_12   |Max. monthly air temperature for December         |average          |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_prec_sum.nc               |1       |prec_sum_01   |Sum of monthly precipitation for January          |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |2       |prec_sum_02   |Sum of monthly precipitation for February         |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |3       |prec_sum_03   |Sum of monthly precipitation for March            |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |4       |prec_sum_04   |Sum of monthly precipitation for April            |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |5       |prec_sum_05   |Sum of monthly precipitation for May              |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |6       |prec_sum_06   |Sum of monthly precipitation for June             |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |7       |prec_sum_07   |Sum of monthly precipitation for July             |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |8       |prec_sum_08   |Sum of monthly precipitation for August           |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |9       |prec_sum_09   |Sum of monthly precipitation for September        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |10      |prec_sum_10   |Sum of monthly precipitation for October          |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |11      |prec_sum_11   |Sum of monthly precipitation for November         |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_sum.nc               |12      |prec_sum_12   |Sum of monthly precipitation for December         |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_tmin_weighted_average.nc  |1       |tmin_wavg_01  |Min. monthly air temperature January              |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |2       |tmin_wavg_02  |Min. monthly air temperature February             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |3       |tmin_wavg_03  |Min. monthly air temperature March                |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |4       |tmin_wavg_04  |Min. monthly air temperature April                |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |5       |tmin_wavg_05  |Min. monthly air temperature May                  |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |6       |tmin_wavg_06  |Min. monthly air temperature June                 |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |7       |tmin_wavg_07  |Min. monthly air temperature July                 |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |8       |tmin_wavg_08  |Min. monthly air temperature August               |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |9       |tmin_wavg_09  |Min. monthly air temperature September            |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |10      |tmin_wavg_10  |Min. monthly air temperature October              |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |11      |tmin_wavg_11  |Min. monthly air temperature November             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmin_weighted_average.nc  |12      |tmin_wavg_12  |Min. monthly air temperature December             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |1       |tmax_wavg_01  |Max. monthly air temperature January              |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |2       |tmax_wavg_02  |Max. monthly air temperature February             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |3       |tmax_wavg_03  |Max. monthly air temperature March                |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |4       |tmax_wavg_04  |Max. monthly air temperature April                |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |5       |tmax_wavg_05  |Max. monthly air temperature May                  |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |6       |tmax_wavg_06  |Max. monthly air temperature June                 |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |7       |tmax_wavg_07  |Max. monthly air temperature July                 |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |8       |tmax_wavg_08  |Max. monthly air temperature August               |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |9       |tmax_wavg_09  |Max. monthly air temperature September            |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |10      |tmax_wavg_10  |Max. monthly air temperature October              |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |11      |tmax_wavg_11  |Max. monthly air temperature November             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_tmax_weighted_average.nc  |12      |tmax_wavg_12  |Max. monthly air temperature December             |weighted average |[°C] * 10      |Int32   |-999   |WorldClim  |
|monthly_prec_weighted_sum.nc      |1       |prec_wsum_01  |Sum of monthly precipitation January              |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |2       |prec_wsum_02  |Sum of monthly precipitation February             |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |3       |prec_wsum_03  |Sum of monthly precipitation March                |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |4       |prec_wsum_04  |Sum of monthly precipitation April                |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |5       |prec_wsum_05  |Sum of monthly precipitation May                  |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |6       |prec_wsum_06  |Sum of monthly precipitation June                 |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |7       |prec_wsum_07  |Sum of monthly precipitation July                 |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |8       |prec_wsum_08  |Sum of monthly precipitation August               |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |9       |prec_wsum_09  |Sum of monthly precipitation September            |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |10      |prec_wsum_10  |Sum of monthly precipitation October              |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |11      |prec_wsum_11  |Sum of monthly precipitation November             |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|monthly_prec_weighted_sum.nc      |12      |prec_wsum_12  |Sum of monthly precipitation December             |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |1       |hydro_avg_01  |Bioclim 1                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |2       |hydro_avg_02  |Bioclim 2                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |3       |hydro_avg_03  |Bioclim 3                                         |average          | * 100         |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |4       |hydro_avg_04  |Bioclim 4                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |5       |hydro_avg_05  |Bioclim 5                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |6       |hydro_avg_06  |Bioclim 6                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |7       |hydro_avg_07  |Bioclim 7                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |8       |hydro_avg_08  |Bioclim 8                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |9       |hydro_avg_09  |Bioclim 9                                         |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |10      |hydro_avg_10  |Bioclim 10                                        |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |11      |hydro_avg_11  |Bioclim 11                                        |average          |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |12      |hydro_avg_12  |Bioclim 12                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |13      |hydro_avg_13  |Bioclim 13                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |14      |hydro_avg_14  |Bioclim 14                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |15      |hydro_avg_15  |Bioclim 15                                        |sum              | * 100         |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |16      |hydro_avg_16  |Bioclim 16                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |17      |hydro_avg_17  |Bioclim 17                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |18      |hydro_avg_18  |Bioclim 18                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_average+sum.nc          |19      |hydro_avg_19  |Bioclim 19                                        |sum              |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |1       |hydro_wavg_01 |Bioclim 1                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |2       |hydro_wavg_02 |Bioclim 2                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |3       |hydro_wavg_03 |Bioclim 3                                         |weighted average | * 100         |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |4       |hydro_wavg_04 |Bioclim 4                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |5       |hydro_wavg_05 |Bioclim 5                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |6       |hydro_wavg_06 |Bioclim 6                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |7       |hydro_wavg_07 |Bioclim 7                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |8       |hydro_wavg_08 |Bioclim 8                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |9       |hydro_wavg_09 |Bioclim 9                                         |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |10      |hydro_wavg_10 |Bioclim 10                                        |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |11      |hydro_wavg_11 |Bioclim 11                                        |weighted average |[°C] * 10      |Int32   |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |12      |hydro_wavg_12 |Bioclim 12                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |13      |hydro_wavg_13 |Bioclim 13                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |14      |hydro_wavg_14 |Bioclim 14                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |15      |hydro_wavg_15 |Bioclim 15                                        |weighted sum     | * 100         |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |16      |hydro_wavg_16 |Bioclim 16                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |17      |hydro_wavg_17 |Bioclim 17                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |18      |hydro_wavg_18 |Bioclim 18                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|hydroclim_weighted_average+sum.nc |19      |hydro_wavg_19 |Bioclim 19                                        |weighted sum     |[mm]           |Float64 |-9999  |WorldClim  |
|landcover_minimum.nc              |1       |lc_min_01     |Evergreen/deciduous needleleaf trees              |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |2       |lc_min_02     |Evergreen broadleaf trees                         |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |3       |lc_min_03     |Deciduous broadleaf trees                         |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |4       |lc_min_04     |Mixed/other trees                                 |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |5       |lc_min_05     |Shrubs                                            |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |6       |lc_min_06     |Herbaceous vegetation                             |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |7       |lc_min_07     |Cultivated and managed vegetation                 |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |8       |lc_min_08     |Regularly flooded shrub/herbaceous vegetation     |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |9       |lc_min_09     |Urban/built-up                                    |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |10      |lc_min_10     |Snow/ice                                          |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |11      |lc_min_11     |Barren lands/sparse vegetation                    |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_minimum.nc              |12      |lc_min_12     |Open water                                        |minimum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |1       |lc_max_01     |Evergreen/deciduous needleleaf trees              |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |2       |lc_max_02     |Evergreen broadleaf trees                         |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |3       |lc_max_03     |Deciduous broadleaf trees                         |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |4       |lc_max_04     |Mixed/other trees                                 |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |5       |lc_max_05     |Shrubs                                            |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |6       |lc_max_06     |Herbaceous vegetation                             |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |7       |lc_max_07     |Cultivated and managed vegetation                 |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |8       |lc_max_08     |Regularly flooded shrub/herbaceous vegetation     |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |9       |lc_max_09     |Urban/built-up                                    |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |10      |lc_max_10     |Snow/ice                                          |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |11      |lc_max_11     |Barren lands/sparse vegetation                    |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_maximum.nc              |12      |lc_max_12     |Open water                                        |maximum          |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |1       |lc_range_01   |Evergreen/deciduous needleleaf trees              |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |2       |lc_range_02   |Evergreen broadleaf trees                         |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |3       |lc_range_03   |Deciduous broadleaf trees                         |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |4       |lc_range_04   |Mixed/other trees                                 |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |5       |lc_range_05   |Shrubs                                            |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |6       |lc_range_06   |Herbaceous vegetation                             |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |7       |lc_range_07   |Cultivated and managed vegetation                 |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |8       |lc_range_08   |Regularly flooded shrub/herbaceous vegetation     |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |9       |lc_range_09   |Urban/built-up                                    |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |10      |lc_range_10   |Snow/ice                                          |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |11      |lc_range_11   |Barren lands/sparse vegetation                    |range            |[%]            |Byte    |-127   |CLC        |
|landcover_range.nc                |12      |lc_range_12   |Open water                                        |range            |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |1       |lc_avg_01     |Evergreen/deciduous needleleaf trees              |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |2       |lc_avg_02     |Evergreen broadleaf trees                         |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |3       |lc_avg_03     |Deciduous broadleaf trees                         |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |4       |lc_avg_04     |Mixed/other trees                                 |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |5       |lc_avg_05     |Shrubs                                            |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |6       |lc_avg_06     |Herbaceous vegetation                             |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |7       |lc_avg_07     |Cultivated and managed vegetation                 |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |8       |lc_avg_08     |Regularly flooded shrub/herbaceous vegetation     |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |9       |lc_avg_09     |Urban/built-up                                    |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |10      |lc_avg_10     |Snow/ice                                          |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |11      |lc_avg_11     |Barren lands/sparse vegetation                    |average          |[%]            |Byte    |-127   |CLC        |
|landcover_average.nc              |12      |lc_avg_12     |Open water                                        |average          |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |1       |lc_wavg_01    |Evergreen/deciduous needleleaf trees              |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |2       |lc_wavg_02    |Evergreen broadleaf trees                         |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |3       |lc_wavg_03    |Deciduous broadleaf trees                         |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |4       |lc_wavg_04    |Mixed/other trees                                 |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |5       |lc_wavg_05    |Shrubs                                            |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |6       |lc_wavg_06    |Herbaceous vegetation                             |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |7       |lc_wavg_07    |Cultivated and managed vegetation                 |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |8       |lc_wavg_08    |Regularly flooded shrub/herbaceous vegetation     |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |9       |lc_wavg_09    |Urban/built-up                                    |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |10      |lc_wavg_10    |Snow/ice                                          |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |11      |lc_wavg_11    |Barren lands/sparse vegetation                    |weighted average |[%]            |Byte    |-127   |CLC        |
|landcover_weighted_average.nc     |12      |lc_wavg_12    |Open water                                        |weighted average |[%]            |Byte    |-127   |CLC        |
|geology_weighted_sum.nc           |1       |geo_wsum_01   |Archean                                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |2       |geo_wsum_02   |Archean, Permian                                  |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |3       |geo_wsum_03   |Cambrian                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |4       |geo_wsum_04   |Cambrian, Precambrian                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |5       |geo_wsum_05   |Cambrian, Proterozoic                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |6       |geo_wsum_06   |Carboniferous                                     |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |7       |geo_wsum_07   |Carboniferous, Devonian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |8       |geo_wsum_08   |Carboniferous, Miocene                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |9       |geo_wsum_09   |Cenozoic                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |10      |geo_wsum_10   |Cenozoic, Mesozoic                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |11      |geo_wsum_11   |Cretaceous                                        |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |12      |geo_wsum_12   |Cretaceous, Carboniferous                         |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |13      |geo_wsum_13   |Cretaceous, Devonian                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |14      |geo_wsum_14   |Cretaceous, Jurassic                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |15      |geo_wsum_15   |Cretaceous, Mississippian                         |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |16      |geo_wsum_16   |Cretaceous, Paleogene, Neogene                    |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |17      |geo_wsum_17   |Cretaceous, Permian                               |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |18      |geo_wsum_18   |Cretaceous, Tertiary                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |19      |geo_wsum_19   |Cretaceous, Triassic                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |20      |geo_wsum_20   |Devonian                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |21      |geo_wsum_21   |Devonian, Cambrian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |22      |geo_wsum_22   |Devonian, Ordovician                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |23      |geo_wsum_23   |Devonian, Proterozoic                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |24      |geo_wsum_24   |Devonian, Silurian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |25      |geo_wsum_25   |Devonian, Silurian, Ordovician                    |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |26      |geo_wsum_26   |Holocene                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |27      |geo_wsum_27   |Ice                                               |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |28      |geo_wsum_28   |Jurassic                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |29      |geo_wsum_29   |Jurassic, Cambrian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |30      |geo_wsum_30   |Jurassic, Carboniferous                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |31      |geo_wsum_31   |Jurassic, Devonian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |32      |geo_wsum_32   |Jurassic, Mississippian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |33      |geo_wsum_33   |Jurassic, Ordovician                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |34      |geo_wsum_34   |Jurassic, Permian                                 |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |35      |geo_wsum_35   |Jurassic, Triassic                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |36      |geo_wsum_36   |Kimberlite                                        |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |37      |geo_wsum_37   |Mesozoic                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |38      |geo_wsum_38   |Mesozoic, Cenozoic                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |39      |geo_wsum_39   |Mesozoic, Paleozoic                               |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |40      |geo_wsum_40   |Mesozoic, Palezoic                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |41      |geo_wsum_41   |Miocene                                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |42      |geo_wsum_42   |Mississippian                                     |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |43      |geo_wsum_43   |Mississippian, Cambrian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |44      |geo_wsum_44   |Mississippian, Devonian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |45      |geo_wsum_45   |Neogene                                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |46      |geo_wsum_46   |Neogene, Paleogene                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |47      |geo_wsum_47   |Ordovician                                        |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |48      |geo_wsum_48   |Ordovician, Cambrian                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |49      |geo_wsum_49   |Ordovician, Proterozoic                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |50      |geo_wsum_50   |Paleogene                                         |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |51      |geo_wsum_51   |Paleogene, Cretaceous                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |52      |geo_wsum_52   |Paleozoic                                         |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |53      |geo_wsum_53   |Paleozoic, Mesozoic                               |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |54      |geo_wsum_54   |Paleozoic, Precambrian                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |55      |geo_wsum_55   |Paleozoic, Proterozoic                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |56      |geo_wsum_56   |Pennsylvanian                                     |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |57      |geo_wsum_57   |Pennsylvanian, Devonian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |58      |geo_wsum_58   |Pennsylvanian, Mississippian                      |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |59      |geo_wsum_59   |Permian                                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |60      |geo_wsum_60   |Permian, Carboniferous                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |61      |geo_wsum_61   |Permian, Devonian                                 |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |62      |geo_wsum_62   |Permian, Mississippian                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |63      |geo_wsum_63   |Permian, Pennsylvanian                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |64      |geo_wsum_64   |Permian, Triassic                                 |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |65      |geo_wsum_65   |Pleistocene                                       |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |66      |geo_wsum_66   |Pliocene, Quaternary                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |67      |geo_wsum_67   |Precambrian                                       |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |68      |geo_wsum_68   |Precambrian, Devonian                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |69      |geo_wsum_69   |Precambrian, Paleozoic                            |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |70      |geo_wsum_70   |Proterozoic                                       |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |71      |geo_wsum_71   |Proterozoic, Archean                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |72      |geo_wsum_72   |Quaternary                                        |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |73      |geo_wsum_73   |Quaternary, Neogne                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |74      |geo_wsum_74   |Quaternary, Tertiary                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |75      |geo_wsum_75   |Salt                                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |76      |geo_wsum_76   |Silurian                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |77      |geo_wsum_77   |Silurian, Cambrian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |78      |geo_wsum_78   |Silurian, Ordovician                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |79      |geo_wsum_79   |Silurian, Proterozoic                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |80      |geo_wsum_80   |Tertiary                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |81      |geo_wsum_81   |Tertiary, Cretaceous                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |82      |geo_wsum_82   |Triassic                                          |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |83      |geo_wsum_83   |Triassic, Carboniferous                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |84      |geo_wsum_84   |Triassic, Devonian                                |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |85      |geo_wsum_85   |Triassic, Mississippian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |86      |geo_wsum_86   |Triassic, Ordovician                              |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |87      |geo_wsum_87   |Triassic, Paleozoic                               |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |88      |geo_wsum_88   |Triassic, Pennsylvanian                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |89      |geo_wsum_89   |Triassic, Permian                                 |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |90      |geo_wsum_90   |Triassic, Proterozoic                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |91      |geo_wsum_91   |Unknown                                           |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|geology_weighted_sum.nc           |92      |geo_wsum_92   |Water                                             |weighted sum     |weighted count |Int32   |-9999  |USGS       |
|soil_minimum.nc                   |1       |soil_min_01   |Soil organic carbon                               |minimum          |[g/kg]         |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |2       |soil_min_02   |Soil pH in H2O                                    |minimum          |pH * 10        |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |3       |soil_min_03   |Sand content mass fraction                        |minimum          |[%]            |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |4       |soil_min_04   |Silt content mass fraction                        |minimum          |[%]            |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |5       |soil_min_05   |Clay content mass fraction                        |minimum          |[%]            |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |6       |soil_min_06   |Coarse fragments (> 2 mm fraction) volumetric     |minimum          |[%]            |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |7       |soil_min_07   |Cation exchange capacity                          |minimum          |[cmol/kg]      |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |8       |soil_min_08   |Bulk density of the fine earth fraction           |minimum          |[kg / m3]      |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |9       |soil_min_09   |Depth to bedrock (R horizon) up to maximum 240 cm |minimum          |[cm]           |Int32   |-999   |ISRIC      |
|soil_minimum.nc                   |10      |soil_min_10   |Probability of occurence (0-100%) of R horizon    |minimum          |[%]            |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |1       |soil_max_01   |Soil organic carbon                               |maximum          |[g/kg]         |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |2       |soil_max_02   |Soil pH in H2O                                    |maximum          |pH * 10        |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |3       |soil_max_03   |Sand content mass fraction                        |maximum          |[%]            |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |4       |soil_max_04   |Silt content mass fraction                        |maximum          |[%]            |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |5       |soil_max_05   |Clay content mass fraction                        |maximum          |[%]            |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |6       |soil_max_06   |Coarse fragments (> 2 mm fraction) volumetric     |maximum          |[%]            |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |7       |soil_max_07   |Cation exchange capacity                          |maximum          |[cmol/kg]      |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |8       |soil_max_08   |Bulk density of the fine earth fraction           |maximum          |[kg / m3]      |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |9       |soil_max_09   |Depth to bedrock (R horizon) up to maximum 240 cm |maximum          |[cm]           |Int32   |-999   |ISRIC      |
|soil_maximum.nc                   |10      |soil_max_10   |Probability of occurence (0-100%) of R horizon    |maximum          |[%]            |Int32   |-999   |ISRIC      |
|soil_range.nc                     |1       |soil_range_01 |Soil organic carbon                               |range            |[g/kg]         |Int32   |-999   |ISRIC      |
|soil_range.nc                     |2       |soil_range_02 |Soil pH in H2O                                    |range            |pH * 10        |Int32   |-999   |ISRIC      |
|soil_range.nc                     |3       |soil_range_03 |Sand content mass fraction                        |range            |[%]            |Int32   |-999   |ISRIC      |
|soil_range.nc                     |4       |soil_range_04 |Silt content mass fraction                        |range            |[%]            |Int32   |-999   |ISRIC      |
|soil_range.nc                     |5       |soil_range_05 |Clay content mass fraction                        |range            |[%]            |Int32   |-999   |ISRIC      |
|soil_range.nc                     |6       |soil_range_06 |Coarse fragments (> 2 mm fraction) volumetric     |range            |[%]            |Int32   |-999   |ISRIC      |
|soil_range.nc                     |7       |soil_range_07 |Cation exchange capacity                          |range            |[cmol/kg]      |Int32   |-999   |ISRIC      |
|soil_range.nc                     |8       |soil_range_08 |Bulk density of the fine earth fraction           |range            |[kg / m3]      |Int32   |-999   |ISRIC      |
|soil_range.nc                     |9       |soil_range_09 |Depth to bedrock (R horizon) up to maximum 240 cm |range            |[cm]           |Int32   |-999   |ISRIC      |
|soil_range.nc                     |10      |soil_range_10 |Probability of occurence (0-100%) of R horizon    |range            |[%]            |Int32   |-999   |ISRIC      |
|soil_average.nc                   |1       |soil_avg_01   |Soil organic carbon                               |average          |[g/kg]         |Int32   |-999   |ISRIC      |
|soil_average.nc                   |2       |soil_avg_02   |Soil pH in H2O                                    |average          |pH * 10        |Int32   |-999   |ISRIC      |
|soil_average.nc                   |3       |soil_avg_03   |Sand content mass fraction                        |average          |[%]            |Int32   |-999   |ISRIC      |
|soil_average.nc                   |4       |soil_avg_04   |Silt content mass fraction                        |average          |[%]            |Int32   |-999   |ISRIC      |
|soil_average.nc                   |5       |soil_avg_05   |Clay content mass fraction                        |average          |[%]            |Int32   |-999   |ISRIC      |
|soil_average.nc                   |6       |soil_avg_06   |Coarse fragments (> 2 mm fraction) volumetric     |average          |[%]            |Int32   |-999   |ISRIC      |
|soil_average.nc                   |7       |soil_avg_07   |Cation exchange capacity                          |average          |[cmol/kg]      |Int32   |-999   |ISRIC      |
|soil_average.nc                   |8       |soil_avg_08   |Bulk density of the fine earth fraction           |average          |[kg / m3]      |Int32   |-999   |ISRIC      |
|soil_average.nc                   |9       |soil_avg_09   |Depth to bedrock (R horizon) up to maximum 240 cm |average          |[cm]           |Int32   |-999   |ISRIC      |
|soil_average.nc                   |10      |soil_avg_10   |Probability of occurence (0-100%) of R horizon    |average          |[%]            |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |1       |soil_wavg_01  |Soil organic carbon                               |weighted average |[g/kg]         |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |2       |soil_wavg_02  |Soil pH in H2O                                    |weighted average |pH * 10        |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |3       |soil_wavg_03  |Sand content mass fraction                        |weighted average |[%]            |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |4       |soil_wavg_04  |Silt content mass fraction                        |weighted average |[%]            |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |5       |soil_wavg_05  |Clay content mass fraction                        |weighted average |[%]            |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |6       |soil_wavg_06  |Coarse fragments (> 2 mm fraction) volumetric     |weighted average |[%]            |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |7       |soil_wavg_07  |Cation exchange capacity                          |weighted average |[cmol/kg]      |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |8       |soil_wavg_08  |Bulk density of the fine earth fraction           |weighted average |[kg / m3]      |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |9       |soil_wavg_09  |Depth to bedrock (R horizon) up to maximum 240 cm |weighted average |[cm]           |Int32   |-999   |ISRIC      |
|soil_weighted_average.nc          |10      |soil_wavg_10  |Probability of occurence (0-100%) of R horizon    |weighted average |[%]            |Int32   |-999   |ISRIC      |
|quality_control.nc                |1       |missing_cells |Cells filled with max. neighbour value            |-                |-              |Byte    |-127   |-          |
|quality_control.nc                |2       |cells_removed |Cells that were removed manually                  |-                |-              |Byte    |-127   |-          |
|==================================|========|==============|==================================================|=================|===============|========|=======|===========|

