library(terra)
library(sf)
library(tidymodels)

# Reading data ####

frame <- st_read("output/Skrim/modeling.gpkg", layer="dataframe") |> 
  mutate(across(.cols = c(ar5cover, ar5soil, dmkdepth), .fns = as.factor)) |> 
  filter(!(ar5cover %in% c(11,12)))
predictors <- rast("output/Skrim/predictors.tif")

plot(predictors$elevation)
plot(frame[,"depth_cm"],add=T)

predictivedomain <- st_read("data/Skrim/Skrim-site.gpkg", "mask_predictivedomain")
plot(predictivedomain)

# RF of full training data ####

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

fit_remotesensing %>% 
  augment(new_data = frame) %>% 
  probably::cal_plot_regression(depth_cm, .pred, smooth = FALSE)

# Uncertainty with QRF ####

mod_qrf <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", quantreg = TRUE, seed = 345) %>% 
  set_mode("regression")
workflow_remotesensing_uncertainty <- 
  workflow() %>% 
  add_model(mod_qrf) %>% 
  add_recipe(recipe_remotesensing)

set.seed(345)
fit_remotesensing_uncertainty <- 
  workflow_remotesensing_uncertainty %>% 
  fit(frame)

# Note: no tidymodels-native syntax for quantile predictions yet (https://github.com/tidymodels/parsnip/issues/119)
# Workaround derived from https://github.com/tidymodels/probably/issues/131#issue-2137662307
quant_predict <- function(fit, new_data, level = 0.9) {
  alpha <- (1 - level)
  quant_pred <- predict(fit, new_data, type = "quantiles", 
                        quantiles = c(alpha / 2, 1 - (alpha / 2)))
  quant_pred <- dplyr::as_tibble(quant_pred)
  quant_pred <- stats::setNames(quant_pred, c(".pred_lower", ".pred_upper"))
  quant_pred
}

frame_baked <- workflows::extract_recipe(fit_remotesensing_uncertainty) %>% 
  bake(frame)
quantiles_remotesensing <- quant_predict(fit_remotesensing_uncertainty$fit$fit$fit, 
                                         frame_baked) %>% 
  mutate(observed = frame$depth_cm, .before = 1)

quantiles_remotesensing %>% 
  arrange(observed) %>% 
  rowid_to_column() %>%
  select(rowid, .pred_lower, .pred_upper, observed) %>%
  pivot_longer(cols = !matches("rowid"), names_to = "what", values_to = "value") %>% 
  ggplot(aes(x= rowid, y = value, color = what)) +
  geom_point()
quantiles_remotesensing <- quantiles_remotesensing %>% 
  mutate(coverage = ifelse(observed >= .pred_lower & 
                             observed <= .pred_upper, 1, 0))
quantiles_remotesensing %>% 
  pull(coverage) %>%
  mean()

# Variable importance ####

fit_remotesensing
import_perm_ranger <- fit_remotesensing %>% 
  extract_fit_parsnip() %>% 
  vip::vi(scale = TRUE)
vip::vip(import_perm_ranger, num_features = 25)

pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, new_data = newdata)$.pred  
}
# pfun(fit_remotesensing, frame)
ftnames <- fit_remotesensing |> 
  extract_recipe() |> 
  summary() |> 
  filter(role == "predictor") |> 
  pull(variable)

set.seed(403)  
import_perm_vip <- fit_remotesensing %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "permute", feature_names = ftnames, 
          train = st_drop_geometry(frame), 
          target = "depth_cm", metric = "rmse",
          pred_wrapper = pfun, nsim = 10, scale = TRUE)
vip::vip(import_perm_vip, geom = "boxplot")
vip::vip(import_perm_vip, num_features = 25)

import_firm <- fit_remotesensing %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "firm", ice = TRUE, feature_names = ftnames, 
          train = st_drop_geometry(frame), scale = TRUE)
vip::vip(import_firm, num_features = 25)

import_shap <- fit_remotesensing %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "shap", 
          pred_wrapper = pfun,
          feature_names = ftnames, 
          train = st_drop_geometry(frame),
          scale = TRUE)
vip::vip(import_shap, num_features = 25)

bind_rows(perm.ranger = import_perm_ranger, 
          perm.vip = import_perm_vip, 
          firm = import_firm, 
          shap = import_shap,
          .id = 'type') |> 
  readr::write_csv("output/Skrim/variable_importance.csv")

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

# Make CV folds ####

# With kNNDM from CAST into tidymodels workflow
library(CAST) #v.1.0.2

predictorspts <- predictors |> 
  as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[predictivedomain,] |> 
  st_transform(st_crs(frame))
st_crs(frame) == st_crs(predpts)

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

# Error estimation with CV ####

evaluation_metrics <- metric_set(rmse, mae, rsq, ccc)

## Within mires ####

### Remote sensing model ####

workflow_remotesensing |> 
  extract_preprocessor()

set.seed(456)
fit_remotesensing_knndm <- 
  workflow_remotesensing %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_remotesensing_knndm)

preds <- fit_remotesensing_knndm %>% 
  collect_predictions() |> 
  select(id, depth_cm, .pred)
preds %>% 
  probably::cal_plot_regression(depth_cm, .pred, smooth = FALSE)

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
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
# warning: A correlation computation is required, but `estimate` is constant and has 0 standard deviation
# Some folds contain only one dmkdepth level
collect_metrics(fit_dmkintercept_knndm)

# sanity check
fit_dmkintercept <- 
  workflow_dmkintercept %>% 
  fit(frame)
extract_fit_parsnip(fit_dmkintercept)
frame |> 
  group_by(dmkdepth) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) # sanity check

### Terrain model ####

recipe_terrain <- 
  recipe(formula = depth_cm ~ 
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000, 
         data = frame)
prep(recipe_terrain, training = frame)

workflow_terrain <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_terrain)

set.seed(456)
fit_terrain_knndm <- 
  workflow_terrain %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(fit_terrain_knndm)

### Leveraging model ####

recipe_leveraging <- 
  recipe(formula = depth_cm ~ ., data = select(frame, !starts_with(c("source", "ar5")))) |> 
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
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_leveraging_knndm)
collect_predictions(fit_leveraging_knndm) |> 
  drop_na(.pred) |> 
  nrow()

bind_rows(
  remotesensing = collect_metrics(fit_remotesensing_knndm),
  dmkintercept = collect_metrics(fit_dmkintercept_knndm),
  terrain = collect_metrics(fit_terrain_knndm),
  leveraging = collect_metrics(fit_leveraging_knndm), 
  .id = "model") |> 
  readr::write_csv("output/Skrim/modelmetrics.csv")

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

# Uncertainty with CV ####

# Fit the model and get predictions using resampling
quantiles.cv <- map(folds$splits, function(spliti) {
  train_data <- analysis(spliti)
  test_data <- assessment(spliti)
  # Fit QRF model
  set.seed(345)
  qrf_wf <- workflow_remotesensing_uncertainty %>% 
    fit(train_data)
  # Make predictions
  test_data_baked <- workflows::extract_recipe(qrf_wf) %>% 
    bake(test_data)
  predictions <- quant_predict(qrf_wf$fit$fit$fit, test_data_baked)
  # Return results
  bind_cols(predictions, test_data)
}) %>% 
  bind_rows()

# Evaluate prediction intervals
quantiles.cv <- quantiles.cv %>%
  mutate(coverage = ifelse(depth_cm >= .pred_lower & depth_cm <= .pred_upper, 1, 0),
         .before = 1)

mean(quantiles.cv$coverage)
hist(quantiles.cv$.pred_upper - quantiles.cv$.pred_lower)
quantiles.cv %>% 
  arrange(depth_cm) %>% 
  rowid_to_column() %>%
  select(rowid, .pred_lower, .pred_upper, depth_cm) %>%
  pivot_longer(cols = !matches("rowid"), names_to = "what", values_to = "value") %>% 
  ggplot(aes(x= rowid, y = value, color = what)) +
  geom_point()

quantiles.cv <- quantiles.cv %>% 
  st_as_sf(sf_column_name = 'geom')
quantiles.cv %>% 
  filter(coverage == 0) %>% 
  ggplot() +
  geom_sf()

quantiles.cv %>% 
  st_write("output/Skrim/modeling.gpkg", layer = "quantiles.cv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
