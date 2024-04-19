library(tidyverse)
library(terra)
library(sf)
library(mlr3verse)
library(mlr3spatial)

training <- st_read("output/modeling.gpkg", layer="training")
df <- st_drop_geometry(training) |> 
  select(-cell)
plot(training[,"depth_cm"])

# RF with default hyperparameters ####
# Evaluation with random CV

tsk_depth <- as_task_regr(df, target = "depth_cm", id = "depth")

lrn_rf <- lrn("regr.ranger", num.trees = 1000, importance = "permutation")

set.seed(123)
lrn_rf$train(tsk_depth)
lrn_rf$model
plot(lrn_rf$predict(tsk_depth)) + coord_flip()
plot(lrn_rf$predict(tsk_depth), type = "residual")
lrn_rf$oob_error()
sqrt(lrn_rf$oob_error())

predictors <- rast("output/predictors.tif")
prediction <- predict_spatial(predictors, lrn_rf, format = "terra")
plot(prediction)
writeRaster(prediction, "output/RFprediction.tif", overwrite=TRUE)

# Variable importance
library(iml)
depth_x <- tsk_depth$data(cols = tsk_depth$feature_names)
depth_y <- tsk_depth$data(cols = tsk_depth$target_names)
predictor <- Predictor$new(lrn_rf, data = depth_x, y = depth_y)
set.seed(123)
importance = FeatureImp$new(predictor, loss = "mse", n.repetitions = 100)
importance$plot()
tibble(names = names(lrn_rf$model$variable.importance),
       perm.importance = lrn_rf$model$variable.importance) |> 
  arrange(desc(perm.importance))

# 10-repeats 10-fold CV
rcv1010 <- rsmp("repeated_cv", repeats = 10, folds = 10)

set.seed(123)
rr <- resample(tsk_depth, lrn_rf, rcv1010)
rr$aggregate(msr("regr.rmse"))
rr$aggregate(msr("regr.rsq"))


# Tuning ####
# # May be implemented with nested resampling for evaluation. 
# # https://www.tidymodels.org/learn/work/nested-resampling/
# 
# # A mlr task has to be created in order to use the package
# # We make an mlr task with the iris dataset here
# # (Classification task with makeClassifTask, Regression Task with makeRegrTask)
# depthtask.mlr <- mlr::makeRegrTask(data = df, target = "depth_cm")
# depthtask.mlr
# 
# # Tuning process; Tuning measure is mse
# # “final recommended hyperparameter setting is calculated by taking the best 5% of all SMBO iterations” (Probst et al., 2019)
# set.seed(123)
# res <- tuneRanger::tuneRanger(depthtask.mlr, iters = 100, num.trees = 1000)
# 
# # Mean of best 5 % of the results
# res
# 
# # Model with the new tuned hyperparameters
# lrn_rf_tuned <- lrn("regr.ranger", num.trees = 1000, mtry = 15, min.node.size = 2,
#                     sample.fraction = 0.9, importance = "permutation")
# 
# set.seed(123)
# lrn_rf_tuned$train(tsk_depth)
# lrn_rf_tuned$model
# lrn_rf_tuned$oob_error()
# sqrt(lrn_rf_tuned$oob_error())

sessionInfo()