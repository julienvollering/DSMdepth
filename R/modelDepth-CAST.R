library(CAST) #v1.0.2
library(geodata)
library(terra)
library(sf)
library(caret)
library(tmap)


frame <- st_read("output/modeling.gpkg", layer="dataframe")
plot(frame[,"depth_cm"])
df <- st_drop_geometry(frame) |> 
  dplyr::select(!starts_with("source"))
predictors <- rast("output/predictors.tif")

plot(predictors$elevation)
plot(frame[,"depth_cm"],add=T)

# Random CV ####

# note that to use the data for model training we have to get rid of the 
# geometry column of the sf object
set.seed(10) # set seed to reproduce the model
model_default <- caret::train(dplyr::select(df, -depth_cm),
                              dplyr::pull(df, depth_cm),
                              method="rf",tuneGrid=data.frame("mtry"=8),
                              importance=TRUE, ntree=500,
                              trControl=trainControl(method="cv",number=3, savePredictions = "final"))


prediction <- predict(predictors,model_default,na.rm=TRUE)
plot(prediction)
model_default
global_validation(model_default)

# kNNDM CV ####

ar5 <- st_read("data/Orskogfjellet-site.gpkg", layer="fkb_ar5_clipped") 
ar5 <- st_transform(ar5, crs(frame))
ar5.myr <- dplyr::filter(ar5, arealtype == 60) |> 
  st_geometry()
plot(ar5.myr)
sa <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea")
plot(sa)
modeldomain <- st_intersection(sa, ar5.myr) |> 
  st_cast("MULTIPOLYGON") |> 
  st_union() |> 
  st_cast("POLYGON")
plot(modeldomain)
predictorspts <- predictors |> 
  as.points(values = FALSE) |> 
  st_as_sf()
predpts <- predictorspts[modeldomain,] |> 
  st_transform(st_crs(frame))
plot(predpts, pch = '.')
predictors_modeldomain <- mask(predictors, vect(modeldomain))

st_crs(frame) == st_crs(predpts)
crs(frame) == crs(predpts)

set.seed(10)
indices_knndm <- knndm(tpoints = frame, 
                       predpoints = dplyr::slice_sample(predpts, n=1e4), 
                       k=10) 
plot(indices_knndm, type = "simple")

indices_knndm # Important! Varying k gives varying W
str(indices_knndm)
unique(indices_knndm$clusters)

frame_plot <- dplyr::mutate(frame, cvfold = as.factor(indices_knndm$clusters))
ggplot() + geom_sf(data = frame_plot, aes(color=cvfold),size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ ggtitle("spatial fold membership by color")

set.seed(10)
model <- caret::train(dplyr::select(df, -depth_cm),
                      dplyr::pull(df, depth_cm),
                      method="rf",
                      tuneGrid=data.frame("mtry"=8), 
                      importance=TRUE, ntree=500,
                      trControl=trainControl(method="cv",
                                             index = indices_knndm$indx_train,
                                             savePredictions = "final"))
model
global_validation(model)
# MAE: Treats all errors equally, making it less sensitive to outliers.
# RMSE: Squares the errors before averaging them, so it gives more weight to larger errors, making it more sensitive to outliers.

prediction <- predict(predictors_modeldomain,model,na.rm=TRUE)
plot(prediction)

plot(caret::varImp(model))

# Nearest neighbor distances ####

NNDdefault <- geodist(frame,predictors,cvfolds =model_default$control$indexOut)
summary(NNDdefault)
NNDdefault |> 
  dplyr::group_by(what) |> 
  dplyr::summarise(quart1 = quantile(dist, probs = 0.25), 
                   mean = mean(dist), 
                   quart3 = quantile(dist, probs = 0.75))
plot(geodist(frame,predictors,cvfolds =model_default$control$indexOut))+ 
  scale_x_log10(labels=round)
NNDknndm <- geodist(frame, modeldomain = predictors_modeldomain, cvfolds = indices_knndm$indx_test)
summary(NNDknndm)
NNDknndm |> 
  dplyr::group_by(what) |> 
  dplyr::summarise(quart1 = quantile(dist, probs = 0.25), 
                   mean = mean(dist), 
                   quart3 = quantile(dist, probs = 0.75))
plot(geodist(frame, modeldomain = predictors_modeldomain, cvfolds = indices_knndm$indx_test)) + 
  scale_x_log10(labels=round)

# Forward feature selection (skipped) ####

# JV: Implementation of `ffs` does not seem safe with regards to data leakage!
# JV: nested CV is necessary.
# “Inner CV is used to tune models and outer CV is used to determine model performance 
# without bias. Fast filter functions for feature selection are provided and the package 
# ensures that filters are nested within the outer CV loop to avoid information leakage 
# from performance test sets.” (Lewis et al., 2023)
# set.seed(10)
# ffsmodel <- ffs(st_drop_geometry(splotdata)[,predictors],
#                     st_drop_geometry(splotdata)$Species_richness,
#                     method="rf", 
#                     tuneGrid=data.frame("mtry"=2),
#                     verbose=FALSE,
#                     ntree=25, #make it faster for this tutorial
#                     trControl=trainControl(method="cv",
#                                            index = indices_knndm$indx_train,
#                                            savePredictions = "final"))
# ffsmodel
# global_validation(ffsmodel)
# ffsmodel$selectedvars
# 
# plot(ffsmodel)
# 
# prediction_ffs <- predict(predictors_sp, ffsmodel, na.rm=TRUE)
# plot(prediction_ffs)

# Area of Applicability ####

### AOA for which the spatial CV error applies:
AOA <- aoa(predictors_modeldomain, model, LPD = TRUE, verbose=TRUE)

tm_shape(prediction)+
  tm_raster(title="Depth (cm)",style="cont")+
  tm_shape(AOA$AOA)+
  tm_raster(palette=c("1"=NA,"0"="grey"),style="cat",legend.show = FALSE)+
  tm_layout(frame=FALSE,legend.outside = TRUE)+
  tm_add_legend(type="fill",col="grey",border.lwd=0, labels="Outside \nAOA")
    
plot(c(AOA$DI,AOA$LPD))

errormodel_DI <- errorProfiles(model,AOA,variable="DI")
errormodel_LPD <- errorProfiles(model,AOA,variable="LPD")

plot(errormodel_DI) # No discernable relation
plot(errormodel_LPD) # No discernable relation

expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)
plot(expected_error_LPD)

# NNDM CV ####

set.seed(10)
tictoc::tic()
indices_nndm <- nndm(tpoints = frame, 
                     predpoints = dplyr::slice_sample(predpts, n=1e4)) 
tictoc::toc() # Runtime > 20 min
plot(indices_nndm, type = "simple")

