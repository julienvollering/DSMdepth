library(tidyverse)
library(sf)
library(terra)

# Radiometrics #### 
radK10m <- rast("data/RAD_miclev_K_10m.tif") 
radK10m[radK10m < 0] <- NA
plot(radK10m)

radTh10m <- rast("data/RAD_miclev_Th_10m.tif") 
radTh10m[radTh10m < 0] <- NA
plot(radTh10m)

radU10m <- rast("data/RAD_miclev_U_10m.tif") 
radU10m[radU10m < 0] <- NA
plot(radU10m)

radTC10m <- rast("data/RAD_miclev_TC_10m.tif") 
radTC10m[radTC10m < 0] <- NA
plot(radTC10m)

probe <- read_csv("data/depth/concatenatedRawProbe.csv")
probe <- st_as_sf(probe, coords = c('utm32E', 'utm32N'), crs = 25832)
probecells <- extract(radK10m, probe, cells = TRUE, xy=TRUE, ID=FALSE) |> 
  select(cell, x, y) |> 
  distinct(cell, .keep_all = TRUE) |> 
  st_as_sf(coords = c('x','y'), crs=crs(radK10m))

# Simple terrain ####
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

## Depth to water table ####

### Whitebox + GRASS GUI approach ####
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
                         output = "output/DTW/DTM10mD8FAca.tif",
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
# Flow paths needs to be binary (null = no channel, 1 = channel), as start points for the cost function. 
# select channels (above threshold) and transform them into binary variables 
# r.cost can be run with three different methods of identifying the starting point(s)... 
# ...from a raster map. All non-NULL cells are considered to be starting points. 
D8FAca <- rast("output/DTW/DTM10mD8FAca.tif")
FIA <- D8FAca
FIA[FIA < 0.25*1e4] <- NA
FIA[!is.na(FIA)] <- 1

ar1507 <- st_read("data/Basisdata_1507_Alesund_25832_FKB-AR5_FGDB/Basisdata_1507_Alesund_25832_FKB-AR5_FGDB.gdb",
                  layer="fkb_ar5_omrade")
ar1535 <- st_read("data/Basisdata_1535_Vestnes_25832_FKB-AR5_FGDB/Basisdata_1535_Vestnes_25832_FKB-AR5_FGDB.gdb",
                  layer="fkb_ar5_omrade")
ar <- bind_rows(ar1507, ar1535) 
water <- filter(ar, arealtype == 82  | arealtype == 81) |> 
  st_transform(crs = crs(FIA)) |> 
  st_crop(FIA)
FIAwater <- mask(FIA, water, inverse=TRUE, updatevalue=1)
writeRaster(FIAwater, filename = "output/DTW/FIA2500m2.tif", overwrite=TRUE)

DTMepsg32632Slope <- rast("output/DTW/DTMepsg32632Slope.tif")
RAD_miclev_TC_mask <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_miclev_TC_mask")
RAD_miclev_TC_mask_500mbuffer <- st_buffer(RAD_miclev_TC_mask, dist= 500)
DTMepsg32632SlopeMask <- mask(DTMepsg32632Slope, RAD_miclev_TC_mask_500mbuffer)
writeRaster(DTMepsg32632SlopeMask, 
            filename = "output/DTW/DTMepsg32632SlopeMask.tif", overwrite=TRUE)

# GRASS r.cost (via QGIS) crashed with above inputs (FIA2500m2.tif and DTMepsg32632SlopeMask.tif)
# Future workaround?
# Divide study area into overlapping tiles
# Maybe use 'gdistance' in R instead of GRASS r.cost, to script loop in R-environment.


### GRASS from R approach ####
# https://doi.org/10.5281/zenodo.5718133

# # Setup GRASS
# library("rgrass")
# link2GI::linkGRASS(radK10m, ver_select = TRUE, quiet=FALSE)
# link2GI::searchGRASSW(quiet=FALSE)
# 
# write_RAST(dtm1mProjected, vname = "dtm1mProjected")
# # The function ‘r.hydrodem’ removes all depressions (flags = ‘a’) from the DEM 
# # which is necessary for calculating interruption-free flow channels
# execGRASS(
#   cmd ='r.hydrodem', input = 'dtm1mProjected', output = 'filledHydroDem', 
#   flags = c('a','overwrite'))
# # D8 Flow Directions (flags = ‘s’), resulting in flow-direction layer 
# # Additionally, the flow accumulation using the created D8 flow directions is created 
# execGRASS(
#   cmd = 'r.watershed', elevation = 'filledHydroDem', accumulation = 'accum', 
#   flags = c('s', 'a', 'overwrite'))
# # A grid of terrain slope is calculated on the original (not-filled) DEM. 
# # Later, a cost function will be applied to these values, why units are set to [percent]. 
# execGRASS(
#   cmd = 'r.slope.aspect', elevation ='dtm1mProjected', slope = 'slope', format = 'percent', 
#   flags = c('a','overwrite'))
# # A function is defined for efficient calculations: 
# calcDTW <- function(fia) { #fia <- 0.25
#   # For calculating the flow paths, it is necessary to define a threshold (t) 
#   # for the minimal flow initiation area (FIA), meaning how much area needs to 
#   # accumulate downward the slope for resulting in a channel with simulated surface water. 
#   # Commonly, t is set between 0.25 ha and 16 ha (32 ha). 
#   t = fia * 10000/res(dtm1mProjected)[1]^2 
#   # flow accumulation is based on number of cells. 
#   # We need upstream contributing area [m^2]. FIA is corrected by resolution of the DEM. 
#   # Flow paths needs to be binary (0 = no channel, 1 = channel), as start points for the cost function. 
#   # select channels (above threshold) and transform them into binary variables 
#   execGRASS(
#     cmd = 'r.mapcalc', 
#     expression = paste0('flowLines = if(accum >= ', t, ', 1, null())'), 
#     flags = c('overwrite')) 
#   # Finally, DTW calculated as minimum height difference (slope in percent scaled by resolution of the layer) 
#   # between each cell and the flow path (surface water) layer using a cost function. 
#   # The cost function (Awaida and Westervelt, 2020) starts at each point of the plot paths 
#   # and sums up the height difference to each raster cell. 
#   # calculate the least-cost of slope [%], starting from the channels
#   execGRASS(
#     cmd = 'r.cost', 
#     input = 'slope', start_raster = 'flowLines', output = 'cost', null_cost = 0, 
#     flags = c('overwrite','k')) #'k' for with Knight's move for more accurate results
#   
#   # read raster from gdal and save it as object 
#   DTWxha<-rast(read_RAST('cost')) 
#   # since GRASS r.slope.aspect gives a measure in percent, the cost-grid needs to be corrected by resolution and divided by 100 to achieve [m] 
#   DTWxha1<- DTWxha*res(dtm1mProjected)[1]/100 
#   writeRaster(DTWxha1, paste0('output/DTW/DTW_FIA_',fia,'_ha.tif'), overwrite = T)
#   } 
#   # and finally run the function for the the desired FIAs.
#   lapply(c(0.25), calcDTW) 

# Collate all predictors ####

predictors <- c(radK10m, radTh10m, radU10m, radTC10m, 
                elevation, 
                terrainMean1m, 
                terrain10m,
                MRVBF, TWImean5m, TWI10m, TWIbilinear20m, TWIbilinear50m)
names(predictors) <- c('radK', 'radTh', 'radU', 'radTC', 
                       'elevation', 
                       'slope1m', 'TPI1m', 'TRI1m', 'roughness1m',
                       'slope10m', 'TPI10m', 'TRI10m', 'roughness10m',
                       'MRVBF', 'TWI5m', 'TWI10m', 'TWI20m', 'TWI50m')

radextent <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_10m_mask")
predictors.extent <- crop(predictors, radextent)
dtmextent <- ext(dtm1mProjected)
predictors.extent <- crop(predictors, dtmextent)

plot(predictors.extent)
writeRaster(predictors, "output/predictors.tif")
