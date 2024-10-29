library(tidyverse)

# Predictors table ####

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

# sessionInfo ####

sessioninfo::session_info()
