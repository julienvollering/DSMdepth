# Predictors table ####

library(tidyverse)

rast <- terra::rast("output/predictors.tif")
names(rast)

preds.table <- read_csv("ms/tables/predictors.csv")
View(preds.table)

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
