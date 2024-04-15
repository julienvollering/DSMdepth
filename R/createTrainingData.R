library(tidyverse)
library(terra)
library(sf)

# Cleaning probe data ####
probe <- read_csv("data/concatenatedRawProbe.csv")
probe <- st_as_sf(probe, coords = c('utm32E', 'utm32N'), crs = 25832)

filter(probe, is.na(depth_cm)) # 2 of 160 sampling cells have soil altered by construction
filter(probe, !is.na(comment)) |>
  arrange(comment) |> 
  print(n=160)

# Remove depths which have changed since LIDAR and radiometric surveys in 2015. 
# Do not remove depths that are affected by infrastructure already in 2015.
# Checked all points with comments: "infill", "fyllmasse", "gravemasser", "road"
probeclean <- probe |> 
  filter(!is.na(depth_cm)) |> 
  filter(!(from == "dgps2" & ptname == 276 & comment == "fyllmasse, ikke torv")) |> 
  filter(!(from == "dgps2" & ptname == 251 & comment == "road, fyllmasse paa veg"))

# Raster origo and resolution ####
radK10m <- rast("data/RAD_miclev_K_10m.tif") 
plot(radK10m)
radK10m[radK10m < 0] <- NA
plot(radK10m)

# Depth per cell ####

probecells <- extract(radK10m, probeclean, cells = TRUE, xy = TRUE, ID = FALSE) |> 
  select(-RAD_miclev_K_10m)
probes <- bind_cols(select(probeclean, depth_cm), probecells)

# Here, cell depth is simply the mean of measurements within the cell, regardless of their position.
# May be improved to account for measurement locations.
celldepth <- probes |> 
  group_by(cell) |> 
  summarize(depth_cm= mean(depth_cm), cell = first(cell), x=first(x), y=first(y)) |> 
  st_drop_geometry() |> 
  st_as_sf(coords = c('x','y'), crs=st_crs(radK10m))
write_sf(celldepth, "output/modeling.gpkg", layer="celldepth", append = FALSE)

# Join predictors ####
predictors <- rast("output/predictors.tif")
training <- celldepth |> 
  select(-cell) |> 
  bind_cols(extract(predictors, celldepth, ID=FALSE, cell=TRUE)) |> 
  select(depth_cm, cell, geometry, everything())
write_sf(training, "output/modeling.gpkg", layer="training", append = FALSE)
