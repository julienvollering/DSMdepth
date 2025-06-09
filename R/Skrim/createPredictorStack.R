library(tidyverse)
library(sf)
library(terra)

crssite <- st_crs("epsg:25833")
sa <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_cast("POLYGON") |> 
  st_geometry() |> 
  st_transform(crssite)

dtm1m <- rast("data/Skrim/DTMmerged.tif")

# Radiometrics #### 
radK <- rast("data/Skrim/NGU-2013-029/Kong_Rad_Area3_K/Kong_Rad_Area3_K.ERS") 
radTh <- rast("data/Skrim/NGU-2013-029/Kong_Rad_Area3_Th/Kong_Rad_Area3_Th.ERS") 
radU <- rast("data/Skrim/NGU-2013-029/Kong_Rad_Area3_U/Kong_Rad_Area3_U.ERS") 
radTC <- rast("data/Skrim/NGU-2013-029/Kong_Rad_Area3_TC/Kong_Rad_Area3_TC.ERS") 

rad <- c(radK, radTh, radU, radTC) |> 
  crop(y = st_transform(st_buffer(sa, 1000), st_crs(radK)))

# Sensitivity analysis: Compare cubic spline vs bilinear resampling methods
rad10m_cubic <- rad |> 
  project(y = "epsg:25833",  
          method = "cubicspline",
          res = 10,
          origin = origin(dtm1m))

rad10m_bilinear <- rad |> 
  project(y = "epsg:25833",  
          method = "bilinear",
          res = 10,
          origin = origin(dtm1m))

# Use cubic spline as default (original method)
rad10m <- rad10m_cubic
plot(rad10m)

# Correlation analysis between cubic spline and bilinear resampling methods
# Sample points for comparison (using field measurement locations)
probe <- read_csv("data/Skrim/depth_all.csv")
probe <- st_as_sf(probe, coords = c('X', 'Y'), crs = 25833)
probecells <- extract(rad10m_cubic, probe, cells = TRUE, xy=TRUE, ID=FALSE) |> 
  select(cell, x, y) |> 
  distinct(cell, .keep_all = TRUE) |> 
  st_as_sf(coords = c('x','y'), crs=crs(rad10m_cubic))

# Extract values at measurement locations for both methods
rad_cubic_values <- extract(rad10m_cubic, probecells, ID=FALSE)
rad_bilinear_values <- extract(rad10m_bilinear, probecells, ID=FALSE)

# Calculate correlations for each radiometric variable
map2_dfr(rad_cubic_values, rad_bilinear_values, 
         ~tibble(correlation = cor(.x, .y, use="complete.obs")),
         .id = "variable") |> 
  mutate(variable = names(rad10m_cubic))

# Visual comparison
par(mfrow=c(2,2))
walk2(rad_cubic_values, rad_bilinear_values, 
      ~plot(.x, .y, 
            xlab="Cubic spline", ylab="Bilinear",
            main=paste(names(rad_cubic_values)[which(map_lgl(rad_cubic_values, identical, .x))],
                      "\nr =", round(cor(.x, .y, use="complete.obs"), 3))))
par(mfrow=c(1,1))

# Simple terrain ####

# Lidar point density of DTM
sa <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_transform(crssite)
metadata <- st_read("data/Skrim/eksport_926374_20240909/dtm1/metadata/dtm1_Metadata.shp") |> 
  select(Pkttetthet, Aarstall)
st_intersection(sa, metadata) |> 
  mutate(area = units::drop_units(st_area(geom)),
         areaprop = area/sum(area)) |> 
  st_drop_geometry()

# Inspect data
plot(dtm1m)
origin(dtm1m)

## Elevation ####
origin(rad10m)
origin(dtm1m)
elevation <- resample(dtm1m, rad10m, method = "average")

## Mean 1m-slope, -TPI, -TRI, -Roughness ####
terrain1m <- terrain(dtm1m, v = c("slope", 'TPI', 'TRI', 'roughness'), 
                     unit = "degrees", neighbors = 8)
terrainMean1m <- resample(terrain1m, rad10m, method = "average")
plot(terrainMean1m)
terrainMean1m

## 10m-slope, -TPI, -TRI, -Roughness ####
terrain10m <- terrain(elevation, v = c("slope", 'TPI', 'TRI', 'roughness'), 
                      unit = "degrees", neighbors = 8)
plot(terrain10m)
terrain10m

# Geomorphometric/hydrological ####
# Use DTM with larger extent than the extent of radiometric data -- to avoid edge effects 

## Multiresolution Index of Valley Bottom Flatness (MRVBF) ####
# https://saga-gis.sourceforge.io/saga_tool_doc/2.2.0/ta_morphometry_8.html

# SAGA 9.3.1
# [2024-09-10/07:27:01] [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Execution started...
# __________
# [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Parameters:
#   Grid System: 1; 30010x 15010y; 185425.5x 6605995.5y
# Elevation: DTMmerged
# MRVBF: MRVBF
# MRRTF: MRRTF
# Initial Threshold for Slope: 16
# Threshold for Elevation Percentile (Lowness): 0.4
# Threshold for Elevation Percentile (Upness): 0.35
# Shape Parameter for Slope: 4
# Shape Parameter for Elevation Percentile: 3
# Update Views: true
# Classify: false
# Maximum Resolution (Percentage): 100
# 
# step: 1, resolution: 1.00, threshold slope 16.00
# step: 2, resolution: 1.00, threshold slope 8.00
# step: 3, resolution: 3.00, threshold slope 4.00
# step: 4, resolution: 9.00, threshold slope 2.00
# step: 5, resolution: 27.00, threshold slope 1.00
# step: 6, resolution: 81.00, threshold slope 0.50
# step: 7, resolution: 243.00, threshold slope 0.25
# step: 8, resolution: 729.00, threshold slope 0.12
# step: 9, resolution: 2187.00, threshold slope 0.06
# step: 10, resolution: 6561.00, threshold slope 0.03
# step: 11, resolution: 19683.00, threshold slope 0.02
# step: 12, resolution: 59049.00, threshold slope 0.01
# __________
# total execution time: 8196000 milliseconds (02h 16m 36s)
# 
# [2024-09-10/09:43:38] [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Execution succeeded (02h 16m 36s)
MRVBF1m <- rast("output/Skrim/MRVBF1m.tif")
MRVBF <- resample(MRVBF1m, rad10m, method = "average")

## TWI ####
# Derived from minimum 5m resolution

library(whitebox)
whitebox::wbt_init()

# 5 m
template5m <- rast(crs = crs(rad10m), extent=ext(rad10m), resolution=5)
dtm5m <- resample(dtm1m, template5m, method="average")
writeRaster(dtm5m, "output/Skrim/whitebox/DTM5m.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/Skrim/whitebox/DTM5m.tif",
  output = "output/Skrim/whitebox/DTM5mfilled.tif") 
# Elapsed Time (excluding I/O): 7s

wbt_d_inf_flow_accumulation(input = "output/Skrim/whitebox/DTM5mfilled.tif",
                            output = "output/Skrim/whitebox/DTM5mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")
# Elapsed Time (excluding I/O): 3s

wbt_slope(dem = "output/Skrim/whitebox/DTM5mfilled.tif",
          output = "output/Skrim/whitebox/DTM5mSlope.tif",
          units = "degrees")
# Elapsed Time (excluding I/O): 3s

wbt_wetness_index(sca = "output/Skrim/whitebox/DTM5mDinfFAsca.tif",
                  slope = "output/Skrim/whitebox/DTM5mSlope.tif",
                  output = "output/Skrim/whitebox/TWI5m.tif")
# Elapsed Time (excluding I/O): 0.417s

TWI5m <- rast("output/Skrim/whitebox/TWI5m.tif")
TWImean5m <- resample(TWI5m, rad10m, method="average")

# 10 m
writeRaster(elevation, "output/Skrim/whitebox/DTM10m.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/Skrim/whitebox/DTM10m.tif",
  output = "output/Skrim/whitebox/DTM10mfilled.tif") 

wbt_d_inf_flow_accumulation(input = "output/Skrim/whitebox/DTM10mfilled.tif",
                            output = "output/Skrim/whitebox/DTM10mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/Skrim/whitebox/DTM10mfilled.tif",
          output = "output/Skrim/whitebox/DTM10mSlope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/Skrim/whitebox/DTM10mDinfFAsca.tif",
                  slope = "output/Skrim/whitebox/DTM10mSlope.tif",
                  output = "output/Skrim/whitebox/TWI10m.tif")

TWI10m <- rast("output/Skrim/whitebox/TWI10m.tif")

# 20 m
template20m <- rast(crs = crs(rad10m), extent=ext(rad10m), resolution=20)
dtm20m <- resample(dtm1m, template20m, method="average")
writeRaster(dtm20m, "output/Skrim/whitebox/DTM20m.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/Skrim/whitebox/DTM20m.tif",
  output = "output/Skrim/whitebox/DTM20mfilled.tif") 

wbt_d_inf_flow_accumulation(input = "output/Skrim/whitebox/DTM20mfilled.tif",
                            output = "output/Skrim/whitebox/DTM20mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/Skrim/whitebox/DTM20mfilled.tif",
          output = "output/Skrim/whitebox/DTM20mSlope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/Skrim/whitebox/DTM20mDinfFAsca.tif",
                  slope = "output/Skrim/whitebox/DTM20mSlope.tif",
                  output = "output/Skrim/whitebox/TWI20m.tif")

TWI20m <- rast("output/Skrim/whitebox/TWI20m.tif")
TWIbilinear20m <- resample(TWI20m, rad10m, method="bilinear")

# 50 m
template50m <- rast(crs = crs(rad10m), extent=ext(rad10m), resolution=50)
dtm50m <- resample(dtm1m, template50m, method="average")
writeRaster(dtm50m, "output/Skrim/whitebox/DTM50m.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/Skrim/whitebox/DTM50m.tif",
  output = "output/Skrim/whitebox/DTM50mfilled.tif") 

wbt_d_inf_flow_accumulation(input = "output/Skrim/whitebox/DTM50mfilled.tif",
                            output = "output/Skrim/whitebox/DTM50mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/Skrim/whitebox/DTM50mfilled.tif",
          output = "output/Skrim/whitebox/DTM50mSlope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/Skrim/whitebox/DTM50mDinfFAsca.tif",
                  slope = "output/Skrim/whitebox/DTM50mSlope.tif",
                  output = "output/Skrim/whitebox/TWI50m.tif")

TWI50m <- rast("output/Skrim/whitebox/TWI50m.tif")
TWIbilinear50m <- resample(TWI50m, rad10m, method="bilinear")

## DTW: Depth to water table ####

### Whitebox + ArcGIS Pro GUI approach ####
library(whitebox)
whitebox::wbt_init()
writeRaster(dtm1m, "output/Skrim/DTW/DTM.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/Skrim/DTW/DTM.tif",
  output = "output/Skrim/DTW/DTMfilled.tif") 
# Elapsed Time (excluding I/O): 5min 32.530s

#  https://www.whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#D8FlowAccumulation
wbt_d8_flow_accumulation(input = "output/Skrim/DTW/DTMfilled.tif",
                         output = "output/Skrim/DTW/DTMD8FAca.tif",
                         out_type = "catchment area")
# Elapsed Time (excluding I/O): 27.934s

wbt_slope(dem = "output/Skrim/DTW/DTM.tif",
          output = "output/Skrim/DTW/DTMSlope.tif",
          units = "percent")
# Elapsed Time (excluding I/O): 28.35s

D8FAca <- rast("output/Skrim/DTW/DTMD8FAca.tif")
FIA <- D8FAca

ar5 <- st_read("data/Skrim/Basisdata_3303_Kongsberg_25832_FKB-AR5_FGDB.gdb", 
               layer="fkb_ar5_omrade")
water <- filter(ar5, arealtype == 82  | arealtype == 81) |> 
  st_transform(crs = crs(FIA)) |> 
  st_crop(FIA)

FIAthresholds <- c(0.25, 0.5, 1, 2, 4, 8, 16)*1e4
FIA <- rep(FIA, length(FIAthresholds))
names(FIA) <- paste0("FIA",FIAthresholds,"m2")
for (i in seq_along(FIAthresholds)) {
  x <- FIAthresholds[i]
  FIA[[i]] <- classify(FIA[[i]], matrix(c(0, x, NA,
                                          x, Inf, 1), ncol=3, byrow = TRUE))
}
waterr <- rasterize(water, FIA) # Mask directly with 'water' not working -- terra issue?
FIAwater <- mask(FIA, waterr, inverse=TRUE, updatevalue=1)
plot(FIAwater)
writeRaster(FIAwater, paste0("output/Skrim/DTW/FIA",FIAthresholds,"m2.tif"), overwrite=TRUE)

DTMSlope <- rast("output/Skrim/DTW/DTMSlope.tif")
sa_1kmbuffer <- st_buffer(sa, dist= 1000)
DTMSlopeMask <- mask(DTMSlope, vect(sa_1kmbuffer))
DTMSlopeMaskUnitless <- DTMSlopeMask/100 # Conversion from %
writeRaster(DTMSlopeMaskUnitless, 
            filename = "output/Skrim/DTW/DTMSlopeUnitless.tif", overwrite=TRUE)

# ArcGIS Pro Distance Accumulation (v.3.1.0)

# Input raster or feature source data     \FIA2500m2.tif
# Output distance accumulation raster     \DTW-FIA2500m2.tif
# Input barrier raster or feature data     
# Input surface raster     
# Input cost raster     \DTMSlopeUnitless.tif
# Vertical factor     BINARY 1 -30 30
# Horizontal factor     BINARY 1 45
# Out back direction raster     \DTW-FIA2500m2-dir.tif
# Distance Method     PLANAR

DTWfiles <- c(
  "output/Skrim/DTW/DTW-FIA2500m2.tif",
  "output/Skrim/DTW/DTW-FIA5000m2.tif",
  "output/Skrim/DTW/DTW-FIA10000m2.tif",
  "output/Skrim/DTW/DTW-FIA20000m2.tif",
  "output/Skrim/DTW/DTW-FIA40000m2.tif",
  "output/Skrim/DTW/DTW-FIA80000m2.tif",
  "output/Skrim/DTW/DTW-FIA160000m2.tif"
)
DTW <- rast(DTWfiles)
plot(DTW, range = c(0,10), xlim = c(2e5,2e5+1000), ylim=c(6615e3,6615e3+1000))
DTWmean1m <- resample(DTW, rad10m, method = "average")

# Collate all predictors ####

predictors <- c(rad10m, 
                elevation, 
                terrainMean1m, 
                terrain10m,
                MRVBF, TWImean5m, TWI10m, TWIbilinear20m, TWIbilinear50m, 
                DTWmean1m)
names(predictors) <- c('radK', 'radTh', 'radU', 'radTC', 
                       'elevation', 
                       'slope1m', 'TPI1m', 'TRI1m', 'roughness1m',
                       'slope10m', 'TPI10m', 'TRI10m', 'roughness10m',
                       'MRVBF', 'TWI5m', 'TWI10m', 'TWI20m', 'TWI50m',
                       'DTW2500', 'DTW5000','DTW10000','DTW20000','DTW40000',
                       'DTW80000','DTW160000')

predictors.extent <- crop(predictors, st_bbox(sa))
predictors.masked <- mask(predictors.extent, vect(st_buffer(sa, 1000)))

plot(predictors.masked)
writeRaster(predictors.masked, "output/Skrim/predictors.tif", overwrite=TRUE)
