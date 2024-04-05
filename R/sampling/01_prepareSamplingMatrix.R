library(terra)
library(sf)
library(dplyr)

rad.files <- list.files("data/", pattern = "RAD_miclev_.*10m.tif$",
                        full.names = TRUE)
rad10m <- rast(rad.files)
rad10m.mask <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_10m_mask")
rad10m <- mask(rad10m, rad10m.mask)
plot(rad10m)

ar5 <- st_read("data/Orskogfjellet-site.gpkg", layer="fkb_ar5_clipped")
ar5 <- st_transform(ar5, crs(rad10m))
ar5.myr <- filter(ar5, arealtype == 60)
plot(rad10m["RAD_miclev_TC_10m"])
plot(st_geometry(ar5.myr), add=TRUE)

slope1m <- rast("data/slope.tif")
slope1m

slope10m <- project(slope1m, rad10m, method = "average")
plot(slope10m)

covariates <- c(rad10m, slope10m)
names(covariates) <- gsub("_miclev", "", names(covariates))
names(covariates) <- gsub("_10m", "", names(covariates))
writeRaster(covariates, "output/samplingCovariates.tif")

cost <- rast("data/rwalkcost.tif")
plot(log10(cost))

cost <- project(cost, covariates)
stack <- c(covariates, cost)

sa <- mask(stack, ar5.myr, touches=FALSE)
sa <- trim(sa)
plot(sa)
df <- as.data.frame(sa, xy = TRUE, na.rm=TRUE)
readr::write_csv(df, "output/preparedSamplingMatrix.csv")
