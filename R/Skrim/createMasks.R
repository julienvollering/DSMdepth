library(tidyverse)
library(terra)
library(sf)

crssite <- st_crs("epsg:25833")

# Study area ####

sa <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_geometry() |> 
  st_transform(crssite)
plot(sa)
st_area(sa) |> 
  units::set_units("km^2")

# AR5 myr ####

ar5 <- st_read("data/Skrim/Basisdata_3303_Kongsberg_25832_FKB-AR5_FGDB.gdb", 
               layer="fkb_ar5_omrade")
ar5 <- st_transform(ar5, crssite)
st_intersection(ar5, sa) |> 
  mutate(area = st_area(SHAPE)) |> 
  st_drop_geometry() |> 
  group_by(arealtype) |> 
  summarize(area = units::set_units(sum(area), "km^2")) |> 
  mutate(percent = area / sum((area)) * 100)

ar5.myr <- dplyr::filter(ar5, arealtype == 60) |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_geometry()
st_write(ar5.myr, "data/Skrim/Skrim-site.gpkg", "mask_ar5myr", 
         append = FALSE)

# Bedrock ####

bedrock <- st_read("data/Skrim/Geologi_33_BerggrunnN250/BerggrunnN250.gdb", 
                   layer="BergartFlate_N250") |> 
  select(contains("bergart"))
bedrock <- st_transform(bedrock, crssite) |> 
  group_by(hovedbergart, hovedbergart_navn) |>
  summarize()
st_intersection(bedrock, sa) |> 
  mutate(area = st_area(SHAPE)) |> 
  st_drop_geometry() |> 
  group_by(hovedbergart, hovedbergart_navn) |> 
  summarize(area = units::set_units(sum(area), "km^2"), .groups = "drop") |> 
  mutate(percent = area / sum((area)) * 100) |> 
  arrange(desc(percent))

# Model predictive domain ####

predictivedomain <- st_intersection(sa, ar5.myr) |> 
  st_cast("MULTIPOLYGON") |> 
  st_union() |> 
  st_cast("POLYGON")
plot(predictivedomain)
st_write(predictivedomain, "data/Skrim/Skrim-site.gpkg", "mask_predictivedomain", 
         append = FALSE)

# sessionInfo ####

sessioninfo::session_info()