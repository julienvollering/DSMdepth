library(tidyverse)
library(terra)
library(sf)

depthpts <- read_csv("data/depth/all.csv") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

# Raster origo and resolution ####
radK10m <- rast("data/RAD_miclev_K_10m.tif") 
plot(radK10m)
radK10m[radK10m < 0] <- NA
plot(radK10m)

# Depth per cell ####

depthpts <- st_transform(depthpts, crs(radK10m))
depthcells <- extract(radK10m, depthpts, cells = TRUE, xy = TRUE, ID = FALSE)
depth <- bind_cols(select(depthpts, source, depth_cm), depthcells) |> 
  filter(!is.na(RAD_miclev_K_10m)) |> 
  select(-RAD_miclev_K_10m)

# Here, cell depth is simply the mean of measurements within the cell, regardless of their position.
# May be improved to account for measurement locations.
celldepth <- depth |> 
  group_by(cell) |> 
  st_drop_geometry() |> 
  summarize(depth_cm= mean(depth_cm), cell = first(cell), x=first(x), y=first(y)) |> 
  st_as_sf(coords = c('x','y'), crs=st_crs(radK10m))

# # Here, cell depth is a *weighted* mean of measurements within the cell, weighted by distance to cell center
# depth |> 
#   nest(pts = c(source,depth_cm,geometry)) |> 
#   st_as_sf(coords=c('x','y'), crs=st_crs(radK10m)) |> 
#   slice_sample(n=2) |> 
#   mutate(
#     meandepth = map_dbl(pts, function(sf) {mean(sf$depth_cm)}),
#     pts = mutate(pts, test = "test")
#   )
#   
# myfunction <- function(pts, geometry, ...) {
#   st_distance(pts[[1]], geometry)
# }
# 
# myfunction(temp$pts, temp$geometry)



# Write cell depths ####

write_sf(celldepth, "output/modeling.gpkg", layer="celldepth", append = FALSE)

# Join predictors ####
predictors <- rast("output/predictors.tif")
training <- celldepth |> 
  select(-cell) |> 
  bind_cols(extract(predictors, celldepth, ID=FALSE, cell=TRUE)) |> 
  select(depth_cm, cell, geometry, everything())
write_sf(training, "output/modeling.gpkg", layer="training", append = FALSE)
