library(tidyverse)
library(sf)
library(tidymodels)

# Reading data ####

frame.depth <- st_read("output/Skrim/modeling.gpkg", layer="dataframe") |> 
  mutate(across(.cols = c(ar5cover, ar5soil, dmkdepth), .fns = as.factor)) |> 
  filter(!(ar5cover %in% c(11,12)))
frame.occ <- st_read("output/Skrim/modeling.gpkg", layer="dataframe.occurrence") |> 
  mutate(across(.cols = c(occurrence, ar5cover, ar5soil), .fns = as.factor)) |> 
  filter(!(ar5cover %in% c(11,12)))

# Evaluation with independent data ####

mod_rf <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", importance = "permutation",
             scale.permutation.importance	= TRUE) %>% 
  set_mode("regression")

## test ####

### preds: radiometric + terrain ####

recipe_RT <- 
  recipe(formula = depth_cm ~ 
           radK + radTh + radU + radTC +
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000,
         data = frame.depth)

workflow_RT <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_RT)

set.seed(456)
fit_RT <- 
  workflow_RT %>% 
  fit(frame.depth)

eval_RT <- fit_RT %>% 
  augment(new_data = frame.occ)

eval_RT %>% 
  roc_curve(truth = occurrence, .pred) %>% 
  autoplot()
eval_RT %>% 
  roc_auc(truth = occurrence, .pred)

### preds: terrain ####

recipe_T <- 
  recipe(formula = depth_cm ~ 
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000, 
         data = frame.depth)

workflow_T <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_T)

set.seed(456)
fit_T <- 
  workflow_T %>% 
  fit(frame.depth)

eval_T <- fit_T %>% 
  augment(new_data = frame.occ)

eval_T %>% 
  roc_curve(truth = occurrence, .pred) %>% 
  autoplot()
eval_T %>% 
  roc_auc(truth = occurrence, .pred)

# sessionInfo ####

sessioninfo::session_info()
