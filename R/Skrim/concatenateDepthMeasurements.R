library(tidyverse)
library(sf)

crssite <- st_crs("epsg:25833")

# Probe data ####

probe <- st_read("data/Skrim/Skrim-site.gpkg", "depth")
filter(probe, !is.na(note))
probe <- pivot_longer(probe, cols = starts_with("probe"), 
                      names_to = "replicate",
                      values_to = "depth_cm")
filter(probe, is.na(depth_cm)) 

# GPR data ####

mire1 <- st_read("data/Skrim/GPRpicks.gpkg", "mire1")
mire2 <- st_read("data/Skrim/GPRpicks.gpkg", "mire2")
mire3 <- st_read("data/Skrim/GPRpicks.gpkg", "mire3")
gpr <- bind_rows(mire1, mire2, mire3) |> 
  mutate(depth_cm = Depth*100) |> 
  select(depth_cm)

# Concatenate sources ####

all <- list(probe = select(probe, depth_cm),
            gpr = select(gpr, depth_cm))
all <- all |> 
  map(\(x) st_transform(x, crs = crssite)) |>
  map(\(x) st_sf(st_set_geometry(x, NULL), geometry = st_geometry(x))) |> 
  bind_rows(.id="source")

all |> 
  st_drop_geometry() |> 
  group_by(source) |> 
  summarize(n = n(), meandepth=mean(depth_cm))

st_write(all, "data/Skrim/depth_all.csv", layer_options = "GEOMETRY=AS_XY", append=FALSE)
