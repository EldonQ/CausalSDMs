library(sf)
library(terra)
china_sf <- st_read("chinashp/china.shp", quiet = TRUE)
r_template <- rast(ext(china_sf), res = 0.1, vals = 0)
crs(r_template) <- "EPSG:4326"
grid_coords <- as.data.frame(r_template, xy = TRUE)
print(nrow(grid_coords))
