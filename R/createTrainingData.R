library(tidyverse)
library(terra)
library(sf)

depthpts <- read_csv("data/depth/all.csv") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

# Raster origo and resolution ####
radK10m <- rast("data/RAD_miclev_K_10m.tif") 
radK10m[radK10m < 0] <- NA
plot(radK10m)

# Depth per cell ####

depthpts <- st_transform(depthpts, crs(radK10m))
depthcells <- extract(radK10m, depthpts, cells = TRUE, xy = TRUE, ID = FALSE)
depth <- bind_cols(select(depthpts, source, depth_cm), depthcells) |> 
  filter(!is.na(RAD_miclev_K_10m)) |> 
  select(-RAD_miclev_K_10m)
cellcenters <- depth |> 
  select(x, y) |> 
  st_drop_geometry() |> 
  st_as_sf(coords = c('x','y'), crs = st_crs(depth))
depth <- depth |> 
  mutate(mFromCellCenter = st_distance(depth, cellcenters, by_element = TRUE),
         weight = units::drop_units(1/(mFromCellCenter)^2),
         weighted_depth_cm = depth_cm*weight) |>
  nest(pts = c(source, depth_cm, mFromCellCenter, weight, weighted_depth_cm, geometry)) |>
  mutate(
    depthMean = map_dbl(pts, function(sf) {mean(sf$depth_cm)}),
    sumweights = map_dbl(pts, function(sf) {sum(sf$weight)}),
    sumweighteddepth = map_dbl(pts, function(sf) {sum(sf$weighted_depth_cm)}),
    depthIDW = sumweighteddepth/sumweights,
    sourceProbe = map_lgl(pts, function(sf) {any(sf$source == "probe")}),
    sourceGPR = map_lgl(pts, function(sf) {any(sf$source == "gpr")}),
    sourceWisen2021 = map_lgl(pts, function(sf) {any(sf$source == "wisen2021")}),
    sourceMyrarkivet = map_lgl(pts, function(sf) {any(sf$source == "myrarkivet")})
  ) |> 
  select(-sumweights, -sumweighteddepth)
# depthMean is simply the overall mean of measurements within the cell, regardless of their position.
# depthIDW is the depth at the cell center from IDW of measurements within the cell

plot(depth$depthMean, depth$depthIDW)
cor(depth$depthMean, depth$depthIDW)
celldepth <- depth |> 
  mutate(depth_cm = round(depthIDW)) |>
  select(-cell, -pts, -depthMean, -depthIDW) |> 
  st_as_sf(coords = c('x','y'), crs = st_crs(radK10m))

# Write cell depths ####

write_sf(celldepth, "output/modeling.gpkg", layer="celldepth", append = FALSE)

# Join predictors ####
predictors <- rast("output/predictors.tif")
training <- celldepth |> 
  bind_cols(extract(predictors, celldepth, ID=FALSE, cell=FALSE)) |> 
  select(depth_cm, starts_with("source"), geometry, everything())
write_sf(training, "output/modeling.gpkg", layer="training", append = FALSE)
