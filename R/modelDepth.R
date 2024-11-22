library(tidyverse)
library(sf)
library(tidymodels)

# Reading data ####

frame <- st_read("output/modeling.gpkg", layer="dataframe") |> 
  mutate(across(.cols = c(ar5cover, ar5soil, dmkdepth), .fns = as.factor)) |> 
  filter(!(ar5cover %in% c(11,12)))
predictors <- terra::rast("output/predictors.tif")

predictivedomain <- st_read("data/Orskogfjellet-site.gpkg", "mask_predictivedomain")

# Make CV folds ####

# With kNNDM from CAST into tidymodels workflow
library(CAST) #v.1.0.2

predictorspts <- predictors |> 
  terra::as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[predictivedomain,] |> 
  st_transform(st_crs(frame))
st_crs(frame) == st_crs(predpts)
nrow(frame)/nrow(predpts)

set.seed(123)
predptssample <- slice_sample(predpts, n=1e3)
# Vary number of folds
knndm_list <- c(20,15,10,5) |>
  purrr::map(\(x) knndm(tpoints = frame, 
                        predpoints = predptssample,
                        k = x))
# Pick number of folds that minimizes W
knndm <- knndm_list[[which.min(map_dbl(knndm_list, \(x) x$W))]]
knndm
plot(knndm, type = "simple")
ggplot() + geom_sf(data = mutate(frame, cvfold = as.factor(knndm$clusters)), 
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

mod_rf <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", importance = "permutation",
             scale.permutation.importance	= TRUE) %>% 
  set_mode("regression")

evaluation_metrics <- metric_set(rmse, mae, rsq, ccc)

## Within mires ####

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
         data = frame)

workflow_RT <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_RT)

set.seed(456)
fit_RT_knndm <- 
  workflow_RT %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_RT_knndm)

### preds: DMK ####

recipe_D <- 
  recipe(formula = depth_cm ~ dmkdepth, data = frame) |> 
  step_unknown(dmkdepth)

workflow_D <- 
  workflow() %>% 
  add_model(linear_reg(mode = "regression", engine = "lm")) %>% 
  add_recipe(recipe_D)

set.seed(456)
fit_D_knndm <- 
  workflow_D %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics)
collect_metrics(fit_D_knndm)

# sanity check
fit_D <- 
  workflow_D %>% 
  fit(frame)
extract_fit_parsnip(fit_D)
frame |> 
  group_by(dmkdepth) |> 
  summarize(n = n(), depth_cm = mean(depth_cm)) # sanity check

### preds: terrain ####

recipe_T <- 
  recipe(formula = depth_cm ~ 
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000, 
         data = frame)

workflow_T <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_T)

set.seed(456)
fit_T_knndm <- 
  workflow_T %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_T_knndm)

### preds: terrain + DMK ####

recipe_TD <- 
  recipe(formula = depth_cm ~ 
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000 +
           dmkdepth, 
         data = frame) %>% 
  step_unknown(dmkdepth) |> 
  step_dummy(dmkdepth)

workflow_TD <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_TD)

set.seed(456)
fit_TD_knndm <- 
  workflow_TD %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_TD_knndm)

### preds: radiometric + terrain + DMK ####

recipe_RTD <- 
  recipe(formula = depth_cm ~ 
           radK + radTh + radU + radTC +
           elevation + 
           slope1m + TPI1m + TRI1m + roughness1m + 
           slope10m + TPI10m + TRI10m + roughness10m + 
           MRVBF + 
           TWI5m + TWI10m + TWI20m + TWI50m + 
           DTW2500 + DTW5000 + DTW10000 + DTW20000 + DTW40000 + DTW80000 + DTW160000 +
           dmkdepth, 
         data = frame) %>% 
  step_unknown(dmkdepth) |> 
  step_dummy(dmkdepth)

workflow_RTD <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_RTD)

set.seed(456)
fit_RTD_knndm <- 
  workflow_RTD %>% 
  fit_resamples(
    resamples = folds,
    metrics = evaluation_metrics,
    control = control_resamples(save_pred = TRUE))
collect_metrics(fit_RTD_knndm)

### collect metrics ####

modelmetrics <- bind_rows(
  RadiometricTerrain = collect_metrics(fit_RT_knndm),
  DMK = collect_metrics(fit_D_knndm),
  Terrain = collect_metrics(fit_T_knndm),
  TerrainDMK = collect_metrics(fit_TD_knndm),
  RadiometricTerrainDMK = collect_metrics(fit_RTD_knndm), 
  .id = "model")
ggplot(modelmetrics) +
  geom_linerange(aes(y = model, 
                     x = mean, 
                     xmin = mean - std_err, 
                     xmax = mean + std_err)) +
  geom_point(aes(y = model, x = mean)) +
  facet_wrap(~.metric, nrow = 1, scales = "free_x", axes ="margins")

readr::write_csv(modelmetrics, "output/modelmetrics.csv")

# Best model in CV: terrain + DMK
preds <- fit_TD_knndm %>% 
  collect_predictions() |> 
  select(id, depth_cm, .pred)
preds %>% 
  probably::cal_plot_regression(depth_cm, .pred, smooth = FALSE)
readr::write_csv(preds, "output/calibrationplot.csv")

## Extrapolating beyond mire ####

frame |> 
  st_drop_geometry() |> 
  group_by(ar5cover == 60) |> 
  summarise(n = n(), mean(depth_cm))

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
custom_splits_extrap <- purrr::discard(custom_splits_extrap, \(x) length(x$out_id) == 0)
  
folds_extrap <- manual_rset(splits = custom_splits_extrap, 
                            ids = names(custom_splits_extrap))
folds_extrap

### preds: radiometric + terrain ####

set.seed(456)
fit_RT_extrap <- 
  workflow_RT %>% 
  fit_resamples(
    resamples = folds_extrap,
    metrics = metric_set(rmse, mae),
    control = control_resamples(save_pred = TRUE))
# Warnings for calculating rsq on assessment folds of 1.
collect_metrics(fit_RT_extrap)

RT_extrap <- collect_predictions(fit_RT_extrap) |> 
  select(.pred, depth_cm)
cor.test(RT_extrap$depth_cm, RT_extrap$.pred)

### preds: terrain ####

set.seed(456)
fit_T_extrap <- 
  workflow_T %>% 
  fit_resamples(
    resamples = folds_extrap,
    metrics = metric_set(rmse, mae),
    control = control_resamples(save_pred = TRUE))
# Warnings for calculating rsq on assessment folds of 1.
collect_metrics(fit_T_extrap)

T_extrap <- collect_predictions(fit_T_extrap) |> 
  select(.pred, depth_cm)
cor.test(T_extrap$depth_cm, T_extrap$.pred)

### preds: none (assume 30cm) ####

fit_null_extrap <- frame |> 
  st_drop_geometry() |> 
  filter(ar5cover != 60) |> 
  mutate(.pred = 30,
         .resid = .pred - depth_cm,
         absoluteError = abs(.resid)) |>
  summarise(mae = mean(absoluteError),
            rmse = sqrt(mean(.resid^2))) |> 
  pivot_longer(cols = everything(), names_to = ".metric", values_to = "mean")
fit_null_extrap

### collect metrics ####

modelmetrics.extrap <- bind_rows(
  RadiometricTerrain = collect_metrics(fit_RT_extrap),
  Terrain = collect_metrics(fit_T_extrap),
  null = fit_null_extrap,
  .id = "model")
readr::write_csv(modelmetrics.extrap, "output/modelmetrics-extrapolation.csv")

# Uncertainty with CV ####

mod_qrf <- 
  rand_forest(mtry = NULL, min_n = 5, trees = 1000) %>% 
  set_engine("ranger", quantreg = TRUE, seed = 345) %>% 
  set_mode("regression")
workflow_TD_uncertainty <- 
  workflow() %>% 
  add_model(mod_qrf) %>% 
  add_recipe(recipe_TD)

# Note: no tidymodels-native syntax for quantile predictions yet (https://github.com/tidymodels/parsnip/issues/119)
# Workaround derived from https://github.com/tidymodels/probably/issues/131#issue-2137662307
quant_predict <- function(fit, new_data, level = 0.9) {
  alpha <- (1 - level)
  quant_pred <- predict(fit, new_data, type = "quantiles", 
                        quantiles = c(alpha / 2, 1 - (alpha / 2)))
  quant_pred <- as_tibble(quant_pred)
  quant_pred <- stats::setNames(quant_pred, c(".pred_lower", ".pred_upper"))
  quant_pred
}
# Fit the model and get predictions using resampling
quantiles.cv <- map(folds$splits, function(spliti) {
  train_data <- analysis(spliti)
  test_data <- assessment(spliti)
  # Fit QRF model
  set.seed(345)
  qrf_wf <- workflow_TD_uncertainty %>% 
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
  st_write("output/modeling.gpkg", layer = "quantiles.cv", append = FALSE)

# Full training data ####

# Retrain on full training data 
set.seed(345)
fit_TD <- 
  workflow_TD %>% 
  fit(frame)

fit_TD %>% 
  augment(new_data = frame) %>% 
  probably::cal_plot_regression(depth_cm, .pred, smooth = FALSE)

# Uncertainty ####

set.seed(345)
fit_TD_uncertainty <- 
  workflow_TD_uncertainty %>% 
  fit(frame)

frame_baked <- workflows::extract_recipe(fit_TD_uncertainty) %>% 
  bake(frame)
quantiles_TD <- quant_predict(fit_TD_uncertainty$fit$fit$fit, 
                              frame_baked) %>% 
  mutate(observed = frame$depth_cm, .before = 1)

quantiles_TD %>% 
  arrange(observed) %>% 
  rowid_to_column() %>%
  select(rowid, .pred_lower, .pred_upper, observed) %>%
  pivot_longer(cols = !matches("rowid"), names_to = "what", values_to = "value") %>% 
  ggplot(aes(x= rowid, y = value, color = what)) +
  geom_point()
quantiles_TD <- quantiles_TD %>% 
  mutate(coverage = ifelse(observed >= .pred_lower & 
                             observed <= .pred_upper, 1, 0))
quantiles_TD %>% 
  pull(coverage) %>%
  mean()

# Feature selection ####

cormat <- bake(prep(recipe_TD, frame), frame) |> 
  select(!depth_cm) |> 
  cor()
corrplot::corrplot(cormat, method = 'number', type = 'upper', diag = TRUE,
                   order = 'hclust', number.cex = 0.5, addCoefasPercent=TRUE)

recipe_TD_uncorr <- recipe_TD %>% 
  step_corr(all_predictors(), threshold = 0.7, method = "pearson") # Of two variables with super-cutoff correlation, that with smallest mean absolute correlation is kept

workflow_TD_uncorr <- 
  workflow() %>% 
  add_model(mod_rf) %>% 
  add_recipe(recipe_TD_uncorr)

set.seed(456)
fit_TD_uncorr <- 
  workflow_TD_uncorr %>% 
  fit(frame)

vars_sel <- extract_fit_engine(fit_TD_uncorr)$forest$independent.variable.names
vars_diff <- setdiff(colnames(cormat), vars_sel)
vars_corr <- map_chr(vars_diff, function(x) {
  correlations <- cormat[x, vars_sel, drop = FALSE]
  colnames(correlations[,which.max(abs(correlations)), drop = FALSE])})
corrKey <- tibble(removed = vars_diff, correlated = vars_corr) |> 
  group_by(correlated) |>
  summarize(removed = paste(removed, collapse = ", "), .groups = 'drop')

## Variable importance ####

pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, new_data = newdata)$.pred  
}
# pfun(fit_TD_uncorr, frame)
ftnames <- fit_TD_uncorr |> 
  extract_recipe() |> 
  summary() |> 
  filter(role == "predictor") |> 
  pull(variable)

set.seed(403)  
import_perm_vip <- fit_TD_uncorr %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "permute", feature_names = ftnames, 
          train = bake(extract_recipe(fit_TD_uncorr), frame), 
          target = "depth_cm", metric = "rmse",
          pred_wrapper = pfun, nsim = 10, scale = TRUE)
vip::vip(import_perm_vip, geom = "boxplot")
vip::vip(import_perm_vip, num_features = 30)

import_firm <- fit_TD_uncorr %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "firm", ice = TRUE, feature_names = ftnames, 
          train = bake(extract_recipe(fit_TD_uncorr), frame), 
          scale = TRUE)
vip::vip(import_firm, num_features = 30)

import_shap <- fit_TD_uncorr %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "shap", 
          pred_wrapper = pfun,
          feature_names = ftnames, 
          train = bake(extract_recipe(fit_TD_uncorr), frame),
          scale = TRUE)
vip::vip(import_shap, num_features = 30)

bind_rows(perm.vip = import_perm_vip, 
          firm = import_firm, 
          shap = import_shap,
          .id = 'type') |> 
  left_join(corrKey, by = c("Variable" = "correlated")) |>
  readr::write_csv("output/variable_importance.csv")

## Partial dependence plots ####
# https://www.tmwr.org/explain#building-global-explanations-from-local-explanations

featuresquant <- import_shap %>% 
  filter(!grepl("dmkdepth", Variable)) |>
  pull(Variable)
featuresqual <- import_shap %>% 
  filter(grepl("dmkdepth", Variable)) |>
  pull(Variable)

# library(DALEXtra)
# explainer_rf <-
#   explain_tidymodels(
#     fit_TD_uncorr,
#     data = as_tibble(frame),
#     y = frame$depth_cm,
#     label = "random forest")
# set.seed(1805)
# pdp <- model_profile(explainer_rf, N = 500, variables = featuresquant)
# plot(pdp)
# 
# library(vivid)
# pdpVars(data = bake(extract_recipe(fit_TD_uncorr), frame),
#         fit = extract_fit_engine(fit_TD_uncorr),
#         response = "depth_cm",
#         vars = c(featuresquant, featuresqual),
#         nIce = 100)

X <- bake(extract_recipe(fit_TD_uncorr), frame) %>% 
  select(-depth_cm)
y <- frame$depth_cm              
predictor <- iml::Predictor$new(
  model = extract_fit_parsnip(fit_TD_uncorr)$fit,
  data = X,
  y = y)
pdpice <- iml::FeatureEffects$new(predictor, 
                             features = c(featuresquant, featuresqual), 
                             method = "pdp+ice",
                             grid.size = 40)

observed <- extract_recipe(fit_TD_uncorr) |> 
  bake(frame) %>% 
  dplyr::select(-depth_cm) %>% 
  pivot_longer(cols = everything(), names_to = "feature", values_to = ".borders") %>% 
  mutate(.type = "observed")
as.list(pdpice$results) %>% 
  bind_rows(.id = "feature") %>% 
  bind_rows(observed) %>%
  readr::write_csv("output/pdpice.csv")

# Residual spatial structure ####

frame.resid <- fit_TD |> 
  augment(new_data = frame) |> 
  dplyr::select(.pred, depth_cm, geom) |> 
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
