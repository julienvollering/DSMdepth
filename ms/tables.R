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
    grepl("djup", dmkdepth, ignore.case = TRUE) ~ "deep (>100 cm)",
    grepl("grunn", dmkdepth, ignore.case = TRUE) ~ "shallow (<100 cm)",
    TRUE ~ "unknown"
  )) |> 
  filter(!(ar5cover %in% c("11", "12"))) |> 
  mutate(ar5cover = case_when(
    ar5cover == "21" ~ "Agricultural",
    ar5cover == "30" ~ "Forest",
    ar5cover == "50" ~ "Open upland",
    ar5cover == "60" ~ "Peatland",
    TRUE ~ "unknown"
  ))
overall <- sites |> 
  group_by(site, TRUE) |> 
  summarize(n = n(), depth_mean = mean(depth_cm), depth_sd = sd(depth_cm)) |> 
  mutate(percent = (n / sum(n))*100) |>
  pivot_wider(names_from = 'site', 
              values_from = c(n, percent, depth_mean, depth_sd), 
              names_vary = "slowest")
byAR5 <- sites |> 
  group_by(site, ar5cover) |> 
  summarize(n = n(), depth_mean = mean(depth_cm), depth_sd = sd(depth_cm)) |> 
  mutate(percent = (n / sum(n))*100) |>
  pivot_wider(names_from = 'site', 
              values_from = c(n, percent, depth_mean, depth_sd), 
              names_vary = "slowest")
byDMK <- sites |> 
  group_by(site, dmkdepth) |> 
  summarize(n = n(), depth_mean = mean(depth_cm), depth_sd = sd(depth_cm)) |> 
  mutate(percent = (n / sum(n))*100) |>
  pivot_wider(names_from = 'site', 
              values_from = c(n, percent, depth_mean, depth_sd), 
              names_vary = "slowest")
bind_rows(overall, byAR5, byDMK) |> 
  mutate(across(where(is.numeric), ~ round(.x, 0))) |>
  mutate(depth_skrim = paste0(depth_mean_skrim, " (", depth_sd_skrim, ")", " cm"),
         depth_orskog = paste0(depth_mean_orskog, " (", depth_sd_orskog, ")", " cm"),
         class = coalesce(ar5cover, dmkdepth),
  ) |>
  select(class, 
         n_skrim, percent_skrim, depth_skrim,
         n_orskog, percent_orskog, depth_orskog) |> 
  mutate(across(starts_with("depth_"), ~ ifelse(grepl("NA", .x), NA, .x))) |> 
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
