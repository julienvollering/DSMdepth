library(terra)
library(sf)
library(tidymodels)

# Reading data ####

frame <- st_read("output/modeling.gpkg", layer="dataframe") |> 
  mutate(across(.cols = c(ar5cover, ar5soil, dmkdepth), .fns = as.factor)) |> 
  filter(!(ar5cover %in% c(11,12)))
predictors <- rast("output/predictors.tif")

plot(predictors$elevation)
plot(frame[,"depth_cm"],add=T)

predictivedomain <- st_read("data/Orskogfjellet-site.gpkg", "mask_predictivedomain")
plot(predictivedomain)

# RF with default hyperparameters ####

glimpse(frame)
frame |> 
  pull(depth_cm) |> 
  summary()

mod_rf <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", importance = "permutation",
             scale.permutation.importance	= TRUE) %>% 
  set_mode("regression")

recipe_remotesensing <- 
  recipe(formula = depth_cm ~ ., data = select(frame, !starts_with("source"))) |> 
  remove_role(ar5cover, old_role = "predictor") |> 
  remove_role(ar5soil, old_role = "predictor") |> 
  remove_role(dmkdepth, old_role = "predictor") |> 
  remove_role(geom, old_role = "predictor")
  
recipe_remotesensing |> 
  summary() |> 
  print(n=30)

workflow_remotesensing <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_remotesensing)

set.seed(345)
fit_remotesensing <- 
  workflow_remotesensing %>% 
  fit(frame)

# Variable importance ####

fit_remotesensing
import_perm <- fit_remotesensing %>% 
  extract_fit_parsnip() %>% 
  vip::vi()
vip::vip(import_perm, num_features = 25)

ftnames <- fit_remotesensing |> 
  extract_recipe() |> 
  summary() |> 
  filter(role == "predictor") |> 
  pull(variable)
import_firm <- fit_remotesensing %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "firm", feature_names = ftnames, 
          train = st_drop_geometry(frame))
vip::vip(import_firm, num_features = 25)

pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, new_data = newdata)$.pred  
}
# pfun(fit_remotesensing, frame)
import_shap <- fit_remotesensing %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "shap", 
          pred_wrapper = pfun,
          feature_names = ftnames, 
          train = st_drop_geometry(frame))
vip::vip(import_shap, num_features = 25)

# Error estimation ####

evaluation_metrics <- metric_set(rmse, mae, rsq)

# With kNNDM from CAST into tidymodels workflow
library(CAST) #v.1.0.2

predictorspts <- predictors |> 
  as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[predictivedomain,] |> 
  st_transform(st_crs(frame))
crs(frame) == crs(predpts)

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

## Overall ####

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

### Remote sensing model ####

workflow_remotesensing |> 
  extract_preprocessor()

set.seed(456)
fit_remotesensing_knndm <- 
  workflow_remotesensing %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(fit_remotesensing_knndm)

### DMK-only model ####

recipe_dmkintercept <- 
  recipe(formula = depth_cm ~ dmkdepth, data = frame) |> 
  step_unknown(dmkdepth)

recipe_dmkintercept |> 
  summary()

workflow_dmkintercept <- 
  workflow() %>% 
  add_model(linear_reg(mode = "regression", engine = "lm")) %>% 
  add_recipe(recipe_dmkintercept)

workflow_dmkintercept |> 
  extract_preprocessor()

set.seed(456)
fit_dmkintercept_knndm <- 
  workflow_dmkintercept %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(fit_dmkintercept_knndm)

# sanity check
fit_dmkintercept <- 
  workflow_dmkintercept %>% 
  fit(frame)
extract_fit_parsnip(fit_dmkintercept)
frame |> 
  group_by(dmkdepth) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) # sanity check

### Leveraging model ####

recipe_leveraging <- 
  recipe(formula = depth_cm ~ ., data = select(frame, !starts_with("source"))) |> 
  remove_role(geom, old_role = "predictor") |> 
  step_unknown(dmkdepth) |> 
  step_dummy(dmkdepth)
recipe_leveraging
prep(recipe_leveraging, training = frame)

workflow_leveraging <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_leveraging)

set.seed(456)
fit_leveraging_knndm <- 
  workflow_leveraging %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(fit_leveraging_knndm)

## Extrapolating beyond mire ####

frame |> 
  st_drop_geometry() |> 
  group_by(ar5cover == 60) |> 
  summarise(mean(depth_cm))

# Creating rsample splits
custom_splits_extrap <- list()
for (k in seq_along(unique(knndm$clusters))) {
  custom_splits_extrap[[k]] <- make_splits(
    list(analysis = intersect(knndm$indx_train[[k]], 
                              which(frame$ar5cover == 60)), 
         assessment = intersect(knndm$indx_test[[k]], 
                                which(frame$ar5cover != 60))), 
    frame)
}
names(custom_splits_extrap) <- paste0("Fold", unique(knndm$clusters))
custom_splits_extrap <- discard(custom_splits_extrap, \(x) length(x$out_id) == 0)
  
folds_extrap <- manual_rset(splits = custom_splits_extrap, 
                            ids = names(custom_splits_extrap))
folds_extrap

### Remote sensing model ####

workflow_remotesensing |> 
  extract_preprocessor()

set.seed(456)
fit_remotesensing_extrap <- 
  workflow_remotesensing %>% 
  fit_resamples(
    resamples = folds_extrap,
    metrics = metric_set(rmse, mae),
    control = control_resamples(save_pred = TRUE))
# Warnings for calculating rsq on assessment folds of 1.
collect_metrics(fit_remotesensing_extrap)

remotesensing_extrap <- collect_predictions(fit_remotesensing_extrap) |> 
  select(.pred, depth_cm)
plot(remotesensing_extrap)
cor.test(remotesensing_extrap$depth_cm, remotesensing_extrap$.pred)

### 30cm model ####

frame |> 
  st_drop_geometry() |> 
  filter(ar5cover != 60) |> 
  mutate(.pred = 30,
         .resid = .pred - depth_cm,
         absoluteError = abs(.resid)) |>
  summarise(mae = mean(absoluteError),
            rmse = sqrt(mean(.resid^2))) 

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
#   prediction_sites = slice_sample(predpts, n=1e4), # Issue to reprex: st_union(predictivedomain) causes no points to be excluded
#   autocorrelation_range = NULL,
#   min_analysis_proportion = 0.5
# )
# tictoc::toc() # 1 sec for 100 data pts, 300 sec for 500 pts
# autoplot(get_rsplit(nndm, 1))

# Residual spatial structure ####

frame.resid <- fit_remotesensing |> 
  augment(new_data = frame) |> 
  select(.pred, depth_cm, geom) |> 
  mutate(.resid = .pred - depth_cm) |> 
  st_as_sf(sf_column_name = 'geom', crs = st_crs(frame))

plot(frame.resid[,'.resid'])

library(gstat)
emp_variog <- variogram(.resid ~ 1, cutoff = 1e4, 
                        data = as(frame.resid, "Spatial"))
print(emp_variog)
plot(emp_variog)

emp_variog <- variogram(.resid ~ 1, cutoff = 1000, 
                        data = as(frame.resid, "Spatial"))
print(emp_variog)
plot(emp_variog)

emp_variog <- variogram(.resid ~ 1, cutoff = 200, 
                        data = as(frame.resid, "Spatial"))
print(emp_variog)
plot(emp_variog)

# sessionInfo ####

sessioninfo::session_info()
