library(terra)
library(sf)
library(tidymodels)

# Reading data ####

frame <- st_read("output/modeling.gpkg", layer="dataframe")
plot(frame[,"depth_cm"])
predictors <- rast("output/predictors.tif")

plot(predictors$elevation)
plot(frame[,"depth_cm"],add=T)

predictivedomain <- st_read("data/Orskogfjellet-site.gpkg", "mask_predictivedomain")
plot(predictivedomain)

dmk <- st_read("data/Orskogfjellet-site.gpkg", "dmkmyr") |> 
  st_transform(st_crs(frame)) |> 
  mutate(depth_class = as.factor(depth_class))
plot(dmk[predictivedomain,"depth_class"])
frame <- st_join(frame, dmk[,"depth_class"])

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

# Variable importance ####

rf_fit
import_perm <- rf_fit %>% 
  extract_fit_parsnip() %>% 
  vip::vi()
vip::vip(import_perm, num_features = 25)

ftnames <- rf_fit |> 
  extract_recipe() |> 
  summary() |> 
  filter(role == "predictor") |> 
  pull(variable)
import_firm <- rf_fit %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "firm", feature_names = ftnames, 
          train = st_drop_geometry(frame))
vip::vip(import_firm, num_features = 25)

pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, new_data = newdata)$.pred  
}
# pfun(rf_fit, frame)
import_shap <- rf_fit %>% 
  extract_fit_parsnip() |> 
  vip::vi(method = "shap", 
          pred_wrapper = pfun,
          feature_names = ftnames, 
          train = st_drop_geometry(frame))
vip::vip(import_shap, num_features = 25)

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
#   prediction_sites = slice_sample(predpts, n=1e4), # Issue to reprex: st_union(predictivedomain) causes no points to be excluded
#   autocorrelation_range = NULL,
#   min_analysis_proportion = 0.5
# )
# tictoc::toc() # 1 sec for 100 data pts, 300 sec for 500 pts
# autoplot(get_rsplit(nndm, 1))

## With kNNDM from CAST into tidymodels workflow ####

library(CAST) #v.1.0.2

predictorspts <- predictors |> 
  as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[predictivedomain,]
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

# Residual spatial structure ####

frame.resid <- rf_fit |> 
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
