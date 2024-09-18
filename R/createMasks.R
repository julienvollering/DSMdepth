library(tidyverse)
library(terra)
library(sf)

# Study area ####

ar5 <- st_read("data/Orskogfjellet-site.gpkg", layer="fkb_ar5_clipped")
RAD_10m_mask <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_10m_mask")
DTM_mask <- st_read("data/Haram Skodje Ørskog Vestnes 2pkt 2015/metadata/Haram Skodje Ørskog Vestnes 2pkt 2015_Klippefil.shp") |> 
  st_transform(st_crs(RAD_10m_mask))
sea_mask <- filter(ar5, arealtype == 82) |> 
  st_union() |> 
  st_transform(crs = st_crs(RAD_10m_mask))
sa <- st_difference(st_intersection(RAD_10m_mask, DTM_mask), sea_mask) |> 
  st_geometry()
plot(sa)
st_area(sa) |> 
  units::set_units("km^2")
st_write(sa, "data/Orskogfjellet-site.gpkg", "mask_studyarea", 
         append = FALSE)

# AR5 myr ####

ar5 <- st_transform(ar5, st_crs(RAD_10m_mask))
st_intersection(ar5, sa) |> 
  mutate(area = st_area(geom)) |> 
  st_drop_geometry() |> 
  group_by(arealtype) |> 
  summarize(area = units::set_units(sum(area), "km^2")) |> 
  mutate(percent = area / sum((area)) * 100)

ar5.myr <- dplyr::filter(ar5, arealtype == 60) |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_geometry()
st_write(ar5.myr, "data/Orskogfjellet-site.gpkg", "mask_ar5myr", 
         append = FALSE)

# Model predictive domain ####

predictivedomain <- st_intersection(sa, ar5.myr) |> 
  st_cast("MULTIPOLYGON") |> 
  st_union() |> 
  st_cast("POLYGON")
plot(predictivedomain)
st_write(predictivedomain, "data/Orskogfjellet-site.gpkg", "mask_predictivedomain", 
         append = FALSE)
