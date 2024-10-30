library(tidyverse)
library(terra)
library(sf)

depthpts <- read_csv("data/Skrim/depth_all.csv") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25833)
occpts <- st_read("data/Skrim/Skrim-site.gpkg", "occurrence")

# Raster origo and resolution ####
radK10m <- rast("output/Skrim/predictors.tif")[[1]] 
plot(radK10m)

# Depth frame ####

## Depth per cell ####

crs(depthpts) == crs(radK10m)
depthcells <- extract(radK10m, depthpts, cells = TRUE, xy = TRUE, ID = FALSE)
depth <- bind_cols(select(depthpts, source, depth_cm), depthcells) |> 
  filter(!is.na(radK)) |> 
  select(-radK)
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

## Write cell depths ####

write_sf(celldepth, "output/Skrim/modeling.gpkg", layer="celldepth", append = FALSE)

## Join predictors ####

predictors <- rast("output/Skrim/predictors.tif")
frame <- celldepth |> 
  bind_cols(extract(predictors, celldepth, ID=FALSE, cell=FALSE)) |> 
  select(depth_cm, starts_with("source"), geometry, everything())

## Join auxiliary attributes ####

ar5 <- st_read("data/Skrim/Basisdata_3303_Kongsberg_25832_FKB-AR5_FGDB.gdb", 
               layer="fkb_ar5_omrade") |> 
  st_transform(st_crs(frame))
frame <- st_join(frame, ar5[,c("arealtype", "grunnforhold")]) |> 
  rename(ar5cover = arealtype, ar5soil = grunnforhold)

dmk <- st_read("data/Skrim/3303_25833_dmk_1397be04-607e-49f3-9174-534948c80d2b_SHAPE/3303_25833_dmk_1397be04-607e-49f3-9174-534948c80d2b_SHAPE.shp") |> 
  st_transform(st_crs(frame))
frame <- st_join(frame, dmk[,"myr"]) |>
  mutate(dmkdepth = case_when(
    myr > 30 ~ "djup",
    myr > 10 ~ "grunn",
    TRUE ~ NA
  )) |> 
  select(-myr)

frame <- frame |> 
  mutate(across(.cols = c(ar5cover, ar5soil, dmkdepth), .fns = as.factor))

## Write data frame ####

write_sf(frame, "output/Skrim/modeling.gpkg", layer="dataframe", append = FALSE)

# Occurrence frame ####

## Occurrence per cell ####

crs(occpts) == crs(radK10m)
occcells <- extract(radK10m, occpts, cells = TRUE, xy = TRUE, ID = FALSE)
occcells %>% 
  pull(cell) %>% 
  duplicated() %>% 
  any()

## Join predictors ####

frame.occ <- occpts |> 
  select(occurrence = peat) %>% 
  bind_cols(extract(predictors, occpts, ID=FALSE, cell=FALSE))

## Join auxiliary attributes ####

frame.occ <- st_join(frame.occ, ar5[,c("arealtype", "grunnforhold")]) |> 
  rename(ar5cover = arealtype, ar5soil = grunnforhold)
frame.occ %>% 
  st_drop_geometry() %>% 
  group_by(ar5cover, ar5soil) %>%
  count()

## Write data frame ####

write_sf(frame.occ, "output/Skrim/modeling.gpkg", layer="dataframe.occurrence", 
         append = FALSE)
