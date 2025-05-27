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

# Compare predictive differences ####

library(tidyverse)

orskog <- read_csv("output/pairwise_comparisons.csv") 
skrim <- read_csv("output/Skrim/pairwise_comparisons.csv") 

orskog |> 
  filter(.metric %in% c("rsq", "mae", "ccc")) |> 
  mutate(
    across(c(estimate,std.error,df,statistic), ~ round(.x, 2)),
    adj.p.value = ifelse(adj.p.value < 0.001, "< 0.001", round(adj.p.value, 3)),
    .metric = case_match(.metric, 
                     "rsq" ~ "R-squared", 
                     "mae" ~ "Mean absolute error", 
                     "ccc" ~ "Concordance correlation"),
    significance = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  rename(metric = .metric, comparison = contrast) |>
  rename_with(str_to_title) |>
  nest_by(Metric) %>%
  pwalk(~ write_csv(.y, paste0("ms/tables/pairwise-orskog-", .x, ".csv")))

skrim |> 
  filter(.metric %in% c("rsq", "mae", "ccc")) |> 
  mutate(
    across(c(estimate,std.error,df,statistic), ~ round(.x, 2)),
    adj.p.value = ifelse(adj.p.value < 0.001, "< 0.001", round(adj.p.value, 3)),
    .metric = case_match(.metric, 
                         "rsq" ~ "R-squared", 
                         "mae" ~ "Mean absolute error", 
                         "ccc" ~ "Concordance correlation"),
    significance = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  rename(metric = .metric, comparison = contrast) |>
  rename_with(str_to_title) |>
  nest_by(Metric) %>%
  pwalk(~ write_csv(.y, paste0("ms/tables/pairwise-skrim-", .x, ".csv")))

list(orskog,skrim) |> 
  set_names(c("orskog", "skrim")) |>
  bind_rows(.id = "site") |> 
  filter(contrast %in% c("Radiometric - Terrain", "RadiometricDMK - TerrainDMK"),
         .metric %in% c("ccc","mae","rsq")) |> 
  mutate(
    significance = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ "")
  ) |> 
  select(site, .metric, contrast, adj.p.value, significance)

# sessionInfo ####

sessioninfo::session_info()
