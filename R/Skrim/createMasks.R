library(tidyverse)
library(terra)
library(sf)

crssite <- st_crs("epsg:25833")

# Study area ####

sa <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_geometry() |> 
  st_transform(crssite)
plot(sa)

# AR5 myr ####

ar5 <- st_read("data/Skrim/Basisdata_3303_Kongsberg_25832_FKB-AR5_FGDB.gdb", 
               layer="fkb_ar5_omrade")
ar5 <- st_transform(ar5, crssite)
ar5.myr <- dplyr::filter(ar5, arealtype == 60) |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_geometry()
st_write(ar5.myr, "data/Skrim/Skrim-site.gpkg", "mask_ar5myr", 
         append = FALSE)

# Model predictive domain ####

predictivedomain <- st_intersection(sa, ar5.myr) |> 
  st_cast("MULTIPOLYGON") |> 
  st_union() |> 
  st_cast("POLYGON")
plot(predictivedomain)
st_write(predictivedomain, "data/Skrim/Skrim-site.gpkg", "mask_predictivedomain", 
         append = FALSE)
