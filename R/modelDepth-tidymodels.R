library(terra)
library(sf)
library(tidymodels)

# Reading data ####

## Data frame and predictor stack ####

frame <- st_read("output/modeling.gpkg", layer="dataframe")
plot(frame[,"depth_cm"])
predictors <- rast("output/predictors.tif")

plot(predictors$elevation)
plot(frame[,"depth_cm"],add=T)

dmk <- st_read("data/Orskogfjellet-site.gpkg", "dmkmyr") |> 
  st_transform(st_crs(frame)) |> 
  mutate(depth_class = as.factor(depth_class))
  
plot(dmk[modeldomain,"depth_class"])
frame <- st_join(frame, dmk[,"depth_class"])

## Predictive domain ####

ar5 <- st_read("data/Orskogfjellet-site.gpkg", layer="fkb_ar5_clipped") 
ar5 <- st_transform(ar5, crs(frame))
ar5.myr <- dplyr::filter(ar5, arealtype == 60) |> 
  st_geometry()
plot(ar5.myr)
sa <- st_read("output/modeling.gpkg", "studyarea_mask") |> 
  st_geometry()
plot(sa)
modeldomain <- st_intersection(sa, ar5.myr) |> 
  st_cast("MULTIPOLYGON") |> 
  st_union() |> 
  st_cast("POLYGON")
predictorspts <- predictors |> 
  as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[modeldomain,] |> 
  st_transform(st_crs(frame))
plot(predpts, pch = '.')
predictors_modeldomain <- mask(predictors, vect(modeldomain))
plot(predictors_modeldomain[[1]])

st_crs(frame) == st_crs(predpts)

# RF with default hyperparameters ####

glimpse(frame)
frame |> 
  pull(depth_cm) |> 
  summary()

#cores <- parallel::detectCores()

rf_mod <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", importance = "permutation",
             scale.permutation.importance	= TRUE) %>% 
  set_mode("regression")

rf_recipe <- 
  recipe(formula = depth_cm ~ ., data = select(frame, !starts_with("source"))) |> 
  remove_role(depth_class, old_role = "predictor") |> 
  remove_role(geom, old_role = "predictor")
  
rf_recipe |> 
  summary() |> 
  print(n=28)

rf_workflow <- 
  workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe)

set.seed(345)
rf_fit <- 
  rf_workflow %>% 
  fit(frame)

## Variable importance ####

rf_fit
rf_fit %>% 
  extract_fit_parsnip() %>% 
  vip::vip(num_features = 25)

ftnames <- rf_fit |> 
  extract_recipe() |> 
  summary() |> 
  filter(role == "predictor") |> 
  pull(variable)

# import_perm <- rf_fit %>% 
#   extract_fit_parsnip() |> 
#   vip::vi(method = "permute", feature_names = ftnames, 
#           train = st_drop_geometry(frame),
#           target = "depth_cm")

import_firm <- rf_fit %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "firm", feature_names = ftnames, 
          train = st_drop_geometry(frame))

# import_shap <- rf_fit %>% 
#   extract_fit_parsnip() |> 
#   vip::vi(method = "shap", feature_names = ftnames, 
#           train = st_drop_geometry(frame))

# Error estimation ####

## With tidymodels-native NNDM ####

# library(spatialsample) #v.0.5.1
# 
# data(ames, package = "modeldata")
# ames_sf <- sf::st_as_sf(ames, coords = c("Longitude", "Latitude"), crs = 4326)
# 
# # Using a small subset of the data, to make the example run faster:
# temp <- ames_sf[1:100, ]
# temp2 <- ames_sf[2001:2100, ]
# plot(st_geometry(temp))
# plot(st_geometry(temp2))
# temp3 <- spatial_nndm_cv(temp, temp2)
# autoplot(get_rsplit(temp3, 1))
# 
# set.seed(123)
# tictoc::tic()
# nndm <- spatial_nndm_cv(
#   frame[1:100,],
#   prediction_sites = slice_sample(predpts, n=1e4), # Issue to reprex: st_union(modeldomain) causes no points to be excluded
#   autocorrelation_range = NULL,
#   min_analysis_proportion = 0.5
# )
# tictoc::toc() # 1 sec for 100 data pts, 300 sec for 500 pts
# autoplot(get_rsplit(nndm, 1))

## With kNNDM from CAST into tidymodels workflow ####

library(CAST) #v.1.0.2

set.seed(123)
predptssample <- dplyr::slice_sample(predpts, n=1e3)
# Vary number of folds
knndm_list <- c(20,15,10,5) |>
  purrr::map(\(x) knndm(tpoints = frame, 
                        predpoints = predptssample,
                        k = x))
# Pick number of folds that minimizes W
knndm <- knndm_list[[which.min(map_dbl(knndm_list, \(x) x$W))]]
knndm
plot(knndm, type = "simple")
ggplot() + geom_sf(data = dplyr::mutate(frame, cvfold = as.factor(knndm$clusters)), 
                   aes(color=cvfold), size=0.5, shape=3)

# Creating rsample splits
custom_splits <- list()
for (k in seq_along(unique(knndm$clusters))) {
  custom_splits[[k]] <- make_splits(list(analysis = knndm$indx_train[[k]], 
                                         assessment = knndm$indx_test[[k]]), 
                                    frame)
}

folds <- manual_rset(splits = custom_splits, 
                     ids = paste0("Fold", unique(knndm$clusters)))
folds

# Evaluate with CV folds
evaluation_metrics <- metric_set(rmse, mae, rsq)

rf_workflow |> 
  extract_preprocessor()

set.seed(456)
rf_fit_rs <- 
  rf_workflow %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(rf_fit_rs)

## Evaluate DMK ####

intercept_mod <- 
  linear_reg(mode = "regression", engine = "lm")

intercept_recipe <- 
  recipe(formula = depth_cm ~ depth_class, data = frame) 

intercept_recipe |> 
  summary()

intercept_workflow <- 
  workflow() %>% 
  add_model(intercept_mod) %>% 
  add_recipe(intercept_recipe)

intercept_workflow |> 
  extract_preprocessor()

set.seed(456)
intercept_fit_rs <- 
  intercept_workflow %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(intercept_fit_rs)

intercept_fit <- 
  intercept_workflow %>% 
  fit(frame)
extract_fit_parsnip(intercept_fit)
frame |> 
  group_by(depth_class) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) #sanity check

# Residual structure ####
