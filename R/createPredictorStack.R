library(tidyverse)
library(sf)
library(terra)

# Radiometrics #### 
# Load original 50m resolution radiometric data for resampling comparison
radK50m <- rast("data/NGU-2015-015/Romsdalsfjorden_RAD_miclev_K.ers") 
radTh50m <- rast("data/NGU-2015-015/Romsdalsfjorden_RAD_miclev_Th.ers") 
radU50m <- rast("data/NGU-2015-015/Romsdalsfjorden_RAD_miclev_U.ers") 
radTC50m <- rast("data/NGU-2015-015/Romsdalsfjorden_RAD_miclev_TC.ers") 

rad50m <- c(radK50m, radTh50m, radU50m, radTC50m)
rad50m[rad50m < 0] <- NA

## Sensitivity analysis: Compare cubic spline vs bilinear resampling methods ####
rad10m_cubic <- rad50m |> 
  project(y = "epsg:25832",  
          method = "cubicspline",
          res = 10)

rad10m_bilinear <- rad50m |> 
  project(y = "epsg:25832",  
          method = "bilinear",
          res = 10)

# Use pre-processed 10m data as default (original method)
radK10m <- rast("data/RAD_miclev_K_10m.tif") 
radK10m[radK10m < 0] <- NA

radTh10m <- rast("data/RAD_miclev_Th_10m.tif") 
radTh10m[radTh10m < 0] <- NA

radU10m <- rast("data/RAD_miclev_U_10m.tif") 
radU10m[radU10m < 0] <- NA

radTC10m <- rast("data/RAD_miclev_TC_10m.tif") 
radTC10m[radTC10m < 0] <- NA

plot(radK10m)

# Correlation analysis between cubic spline and bilinear resampling methods
# Sample points for comparison (using field measurement locations)
probe <- read_csv("data/depth/concatenatedRawProbe.csv")
probe <- st_as_sf(probe, coords = c('utm32E', 'utm32N'), crs = 25832)
probecells <- extract(radK10m, probe, cells = TRUE, xy=TRUE, ID=FALSE) |> 
  select(cell, x, y) |> 
  distinct(cell, .keep_all = TRUE) |> 
  st_as_sf(coords = c('x','y'), crs=crs(radK10m))

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

## Dose rate ####
radDR10m <- 
  13.078 * radK10m +
   5.675 * radU10m +
   2.494 * radTh10m
names(radDR10m) <- "rad_DR_10m"
c(radDR10m, radTC10m, radK10m, radU10m, radTh10m) |> 
  pairs(hist=FALSE, cor=TRUE, maxcells=1e4)
cor(values(radDR10m), values(radTC10m), use = "pairwise.complete.obs")

# Simple terrain ####

# Lidar point density of DTM
saland <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarealand")
laserdensity <- rast("data/Haram Skodje Ørskog Vestnes 2pkt 2015/metadata/Haram Skodje Ørskog Vestnes 2pkt 2015_Punkttetthet_DTM.tif")
extract(laserdensity, saland, ID=FALSE) |> 
  summary()

# Inspect data
dtm1m <- rast("data/DTMmerged.tif")
plot(dtm1m)
origin(dtm1m)

# Reproject then upscale rather than upscale then reproject
tictoc::tic()
dtm1mProjected <- project(dtm1m, crs(radK10m), origin = origin(radK10m),
                          res = 1, method = "bilinear")
tictoc::toc() # 267 sec
plot(dtm1mProjected)
writeRaster(dtm1mProjected, "output/DTMepsg32632.tif")
#dtm1mProjected <- rast("output/DTMepsg32632.tif")

## Elevation ####
origin(radK10m)
origin(dtm1mProjected)
elevation <- resample(dtm1mProjected, radK10m, method = "average")
# elevation2 <- aggregate(dtm1mProjected, fact=10, fun = "mean")
# plot(extract(elevation, probecells)[,2], 
#      extract(elevation2, probecells)[,2])

## Mean 1m-slope, -TPI, -TRI, -Roughness ####
terrain1m <- terrain(dtm1mProjected, v = c("slope", 'TPI', 'TRI', 'roughness'), 
                     unit = "degrees", neighbors = 8)
terrainMean1m <- resample(terrain1m, radK10m, method = "average")
plot(terrainMean1m)
terrainMean1m

## 10m-slope, -TPI, -TRI, -Roughness ####
terrain10m <- terrain(elevation, v = c("slope", 'TPI', 'TRI', 'roughness'), 
                      unit = "degrees", neighbors = 8)
plot(terrain10m)
terrain10m

# Comparison Mean 1m- vs 10m-
op <- par(mfrow=c(2,2))
walk2(as_tibble(values(terrainMean1m)[probecells$cell,]),
      as_tibble(values(terrain10m)[probecells$cell,]),
      \(x, y) plot(x, y))
par(op)

# Difference in slope between sampling design and predictor stack.
# For sampling design, slope was calculated in QGIS at 1m resolution 
# then resampled to radiometric raster by method = "average"
# Rationale for change: terrain predictors should all be derived from the same, resampled DTM.
slopeSampling <- rast("data/slope.tif")
plot(extract(terrainMean1m$slope, probecells)$slope, 
     extract(slopeSampling, probecells)$slope)
plot(extract(terrain10m$slope, probecells)$slope, 
     extract(slopeSampling, probecells)$slope)

# Geomorphometric/hydrological ####
# Use DTM with larger extent than the extent of radiometric data -- to avoid edge effects 

## Multiresolution Index of Valley Bottom Flatness (MRVBF) ####
# https://saga-gis.sourceforge.io/saga_tool_doc/2.2.0/ta_morphometry_8.html

# SAGA 9.3.2
# [2024-04-12/09:30:14] [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Execution started...
# __________
# [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Parameters:
#   Grid System: 1; 24000x 21600y; 380001x 6928200.2y
# Elevation: DTMepsg32632
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
# total execution time: 5357000 milliseconds (01h 29m 17s)
# 
# [2024-04-12/10:59:32] [Multiresolution Index of Valley Bottom Flatness (MRVBF)] Execution succeeded (01h 29m 17s)
MRVBF1m <- rast("output/MRVBF1m.tif")
tictoc::tic()
MRVBFprojected <- project(MRVBF1m, crs(radK10m), origin = origin(radK10m),
                          res = 1, method = "bilinear")
tictoc::toc() # 40 sec
MRVBF <- resample(MRVBFprojected, radK10m, method = "average")

## TWI ####
# Derived from minimum 5m resolution

library(whitebox)
whitebox::wbt_init()

# 5 m
template5m <- rast(crs = crs(radK10m), extent=ext(radK10m), resolution=5)
dtm5mProjected <- resample(dtm1mProjected, template5m, method="average")
writeRaster(dtm5mProjected, "output/whitebox/DTM5mepsg32632.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/whitebox/DTM5mepsg32632.tif",
  output = "output/whitebox/DTM5mepsg32632filled.tif") 
# Elapsed Time (excluding I/O): 7s

wbt_d_inf_flow_accumulation(input = "output/whitebox/DTM5mepsg32632filled.tif",
                            output = "output/whitebox/DTM5mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")
# Elapsed Time (excluding I/O): 3s

wbt_slope(dem = "output/whitebox/DTM5mepsg32632filled.tif",
          output = "output/whitebox/DTM5mepsg32632Slope.tif",
          units = "degrees")
# Elapsed Time (excluding I/O): 3s

wbt_wetness_index(sca = "output/whitebox/DTM5mDinfFAsca.tif",
                  slope = "output/whitebox/DTM5mepsg32632Slope.tif",
                  output = "output/whitebox/TWI5m.tif")
# Elapsed Time (excluding I/O): 0.417s

TWI5m <- rast("output/whitebox/TWI5m.tif")
TWImean5m <- resample(TWI5m, radK10m, method="average")
# TWIbilinear5m <- resample(TWI5m, radK10m, method="bilinear")
# plot(extract(TWImean5m, probecells)$TWI5m,
#      extract(TWIbilinear5m, probecells)$TWI5m)

# 10 m
writeRaster(elevation, "output/whitebox/DTM10mepsg32632.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/whitebox/DTM10mepsg32632.tif",
  output = "output/whitebox/DTM10mepsg32632filled.tif")

wbt_d_inf_flow_accumulation(input = "output/whitebox/DTM10mepsg32632filled.tif",
                            output = "output/whitebox/DTM10mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/whitebox/DTM10mepsg32632filled.tif",
          output = "output/whitebox/DTM10mepsg32632Slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/whitebox/DTM10mDinfFAsca.tif",
                  slope = "output/whitebox/DTM10mepsg32632Slope.tif",
                  output = "output/whitebox/TWI10m.tif")

TWI10m <- rast("output/whitebox/TWI10m.tif")

# 20 m
template20m <- rast(crs = crs(radK10m), extent=ext(radK10m), resolution=20)
dtm20mProjected <- resample(dtm1mProjected, template20m, method="average")
writeRaster(dtm20mProjected, "output/whitebox/DTM20mepsg32632.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/whitebox/DTM20mepsg32632.tif",
  output = "output/whitebox/DTM20mepsg32632filled.tif") 

wbt_d_inf_flow_accumulation(input = "output/whitebox/DTM20mepsg32632filled.tif",
                            output = "output/whitebox/DTM20mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/whitebox/DTM20mepsg32632filled.tif",
          output = "output/whitebox/DTM20mepsg32632Slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/whitebox/DTM20mDinfFAsca.tif",
                  slope = "output/whitebox/DTM20mepsg32632Slope.tif",
                  output = "output/whitebox/TWI20m.tif")

TWI20m <- rast("output/whitebox/TWI20m.tif")
TWIbilinear20m <- resample(TWI20m, radK10m, method="bilinear")

# 50 m
template50m <- rast(crs = crs(radK10m), extent=ext(radK10m), resolution=50)
dtm50mProjected <- resample(dtm1mProjected, template50m, method="average")
writeRaster(dtm50mProjected, "output/whitebox/DTM50mepsg32632.tif", overwrite=TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = "output/whitebox/DTM50mepsg32632.tif",
  output = "output/whitebox/DTM50mepsg32632filled.tif") 

wbt_d_inf_flow_accumulation(input = "output/whitebox/DTM50mepsg32632filled.tif",
                            output = "output/whitebox/DTM50mDinfFAsca.tif",
                            out_type = "Specific Contributing Area")

wbt_slope(dem = "output/whitebox/DTM50mepsg32632filled.tif",
          output = "output/whitebox/DTM50mepsg32632Slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "output/whitebox/DTM50mDinfFAsca.tif",
                  slope = "output/whitebox/DTM50mepsg32632Slope.tif",
                  output = "output/whitebox/TWI50m.tif")

TWI50m <- rast("output/whitebox/TWI50m.tif")
TWIbilinear50m <- resample(TWI50m, radK10m, method="bilinear")

## DTW: Depth to water table ####

### Whitebox + ArcGIS Pro GUI approach ####
library(whitebox)
whitebox::wbt_init()
writeRaster(dtm1mProjected, "output/DTW/DTMepsg32632.tif", overwrite=TRUE)

# Purchased license required
# wbt_depth_to_water(
#   dem = "output/DTW/DTMepsg32632.tif",
#   output = "output/DTW/DTW.tif")

wbt_fill_depressions_wang_and_liu(
  dem = "output/DTW/DTMepsg32632.tif",
  output = "output/DTW/DTMepsg32632filled.tif") 
# Elapsed Time (excluding I/O): 2min 27.571s

#  https://www.whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#D8FlowAccumulation
wbt_d8_flow_accumulation(input = "output/DTW/DTMepsg32632filled.tif",
                         output = "output/DTW/DTMD8FAca.tif",
                         out_type = "catchment area")
# Elapsed Time (excluding I/O): 21.934s

wbt_slope(dem = "output/DTW/DTMepsg32632.tif",
          output = "output/DTW/DTMepsg32632Slope.tif",
          units = "percent")
# Elapsed Time (excluding I/O): 28.35s

# For calculating the flow paths, it is necessary to define a threshold (t) 
# for the minimal flow initiation area (FIA), meaning how much area needs to 
# accumulate downward the slope for resulting in a channel with simulated surface water. 
# Commonly, t is set between 0.25 ha and 16 ha.
# We try 0.25 and 4 ha (@agrenEvaluatingDigitalTerrain2014; @schonauerSpatiotemporalPredictionSoil2021)
# Flow paths needs to be binary (null = no channel, 1 = channel), as start points for the cost function. 
# select channels (above threshold) and transform them into binary variables 
# r.cost can be run with three different methods of identifying the starting point(s)... 
# ...from a raster map. All non-NULL cells are considered to be starting points. 
D8FAca <- rast("output/DTW/DTMD8FAca.tif")
FIA <- D8FAca

ar5 <- st_read("data/Orskogfjellet-site.gpkg", layer="fkb_ar5_clipped")
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
writeRaster(FIAwater, paste0("output/DTW/FIA",FIAthresholds,"m2.tif"), overwrite=TRUE)

DTMepsg32632Slope <- rast("output/DTW/DTMepsg32632Slope.tif")
RAD_miclev_TC_mask <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_miclev_TC_mask")
RAD_miclev_TC_mask_1kmbuffer <- st_buffer(RAD_miclev_TC_mask, dist= 1000)
DTMepsg32632SlopeMask <- mask(DTMepsg32632Slope, RAD_miclev_TC_mask_1kmbuffer)
DTMepsg32632SlopeMaskUnitless <- DTMepsg32632SlopeMask/100 # Conversion from %
writeRaster(DTMepsg32632SlopeMaskUnitless, 
            filename = "output/DTW/DTMepsg32632SlopeUnitless.tif", overwrite=TRUE)

# ArcGIS Pro Distance Accumulation (v.3.1.0)

# Input raster or feature source data     \FIA2500m2.tif
# Output distance accumulation raster     \DTW-FIA2500m2.tif
# Input barrier raster or feature data     
# Input surface raster     
# Input cost raster     \DTMepsg32632SlopeUnitless.tif
# Vertical factor     BINARY 1 -30 30
# Horizontal factor     BINARY 1 45
# Out back direction raster     \DTW-FIA2500m2-dir.tif
# Distance Method     PLANAR

# Warning: Cost cells with a zero or negative value were encountered. Not all source cells may have been processed.
# (Elapsed Time: 3 minutes 52 seconds)

DTWfiles <- c(
  "output/DTW/DTW-FIA2500m2.tif",
  "output/DTW/DTW-FIA5000m2.tif",
  "output/DTW/DTW-FIA10000m2.tif",
  "output/DTW/DTW-FIA20000m2.tif",
  "output/DTW/DTW-FIA40000m2.tif",
  "output/DTW/DTW-FIA80000m2.tif",
  "output/DTW/DTW-FIA160000m2.tif"
)
DTW <- rast(DTWfiles)
plot(DTW, range = c(0,10), xlim = c(392100,393000), ylim=c(6933800,6934300))
DTWmean1m <- resample(DTW, radK10m, method = "average")

# Collate all predictors ####

predictors <- c(radK10m, radTh10m, radU10m, radTC10m, 
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

sa <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea")
predictors.extent <- crop(predictors, st_bbox(sa))
predictors.masked <- mask(predictors.extent, sa)

plot(predictors.masked)
writeRaster(predictors.masked, "output/predictors.tif", overwrite=TRUE)
