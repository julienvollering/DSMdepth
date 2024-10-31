# Predictors table ####

library(tidyverse)

rast <- terra::rast("output/predictors.tif")
names(rast)
tribble(
  ~name, ~description,
  "radK", "Potassium ground concentration",        
  "radTh", "Thorium ground concentration",
  "radU" , "Uranium ground concentration",        
  "radTC", "Total count of gamma radiation",        
  "elevation", "Mean elevation",
  "slope1m", "Mean of 1 m slope",      
  "TPI1m", "Mean of 1 m topographic position index",
  "TRI1m", "Mean of 1 m terrain ruggedness index",        
  "roughness1m", "Mean of 1 m roughness",  
  "slope10m", "10 m slope",    
  "TPI10m", "10 m topographic position index",
  "TRI10m", "10 m terrain ruggedness index",
  "roughness10m" , "10 m roughness",
  "MRVBF", "Multi-resolution valley bottom flatness",
  "TWI5m", "Mean of 5 m topographic wetness index",
  "TWI10m", "10 m topographic wetness index",
  "TWI20m", "Bilinear interpolation of 20 m topographic wetness index",
  "TWI50m", "Bilinear interpolation of 50 m topographic wetness index",
  "DTW2500", "Depth-to-water index, flow initiation area of 0.25 ha",
  "DTW5000", "Depth-to-water index, flow initiation area of 0.5 ha",
  "DTW10000", "Depth-to-water index, flow initiation area of 1 ha",
  "DTW20000", "Depth-to-water index, flow initiation area of 2 ha",
  "DTW40000", "Depth-to-water index, flow initiation area of 4 ha",
  "DTW80000", "Depth-to-water index, flow initiation area of 8 ha",
  "DTW160000", "Depth-to-water index, flow initiation area of 16 ha",
  "DMK", "DMK peat depth class, categorical with 3 levels") |> 
  write_csv("ms/tables/predictors.csv")

# Summary cell depths by attributes ####

library(tidyverse)

skrim <- sf::st_read("output/Skrim/modeling.gpkg", layer="dataframe") |> 
  sf::st_drop_geometry()
orskog <- sf::st_read("output/modeling.gpkg", layer="dataframe") |> 
  sf::st_drop_geometry()
sites <- bind_rows(skrim = skrim, orskog = orskog, .id = "site") |> 
  mutate(dmkdepth = case_when(
    dmkdepth == "djup myr" ~ "djup",
    dmkdepth == "grunn myr" ~ "grunn",
    TRUE ~ dmkdepth
  )) |> 
  filter(!(ar5cover %in% c("11", "12")))
overall <- sites |> 
  group_by(site) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) |> 
  pivot_wider(names_from = 'site', values_from = c(n, depth_cm), 
              names_vary = "slowest")
byAR5 <- sites |> 
  group_by(site, ar5cover) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) |> 
  mutate(prop = n / sum(n)) |>
  pivot_wider(names_from = 'site', values_from = c(n, prop, depth_cm), 
              names_vary = "slowest")
byDMK <- sites |> 
  group_by(site, dmkdepth) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) |> 
  mutate(prop = n / sum(n)) |>
  pivot_wider(names_from = 'site', values_from = c(n, prop, depth_cm), 
              names_vary = "slowest")
bind_rows(overall, byAR5, byDMK) |> 
  select(ar5cover, dmkdepth, n_skrim, prop_skrim, depth_cm_skrim,
         n_orskog, prop_orskog, depth_cm_orskog) |> 
  write_csv("ms/tables/celldepth-attribute.csv")

# sessionInfo ####

sessioninfo::session_info()
