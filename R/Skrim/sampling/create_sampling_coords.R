library(iSDM) # source("R/Skrim/iSDM-eSample.R")
library(raster)
library(sp)
library(sf)
?eSample

sa <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline")
sa <- st_transform(sa, crs=25833)

mire <- st_read("data/Skrim/Skrim-site.gpkg", "myrWMS")

preds <- stack(list.files(path = "data/Skrim/", pattern= "aligned_", full.names = TRUE))
plot(preds[[1]])
plot(sa, add=TRUE)
preds_sa <- mask(preds, sa)
plot(preds_sa)

# Select sample coordinates for depth model
preds_mire <- mask(preds_sa, mire)
plot(preds_mire)
sum(complete.cases(values(preds_mire)))
set.seed(123)
dsample <- eSample(envData = preds_mire, 
                   nExpect = 100, 
                   saveShape = FALSE,
                   nf = 2, 
                   lowerLim = 1e-2, upperLim = 1-1e-2)
dcoords <- st_as_sf(dsample$GeoSamples)
dcoords <- st_set_crs(dcoords, 25833)
#st_write(dcoords, "dcoords.csv")

# Select sample coordinates for occurrence model

## Discard criterion based on maximum slope: threshold value? 
exp(3) # 20.08554,  See Gatis et al. 2019 figure 2b
mireslope <- values(preds_mire[["aligned_Slope"]])
quantile(mireslope, probs = c(0.9,0.99,0.999), na.rm=TRUE)
mireslope_cdf <- ecdf(mireslope)
mireslope_cdf(20)

sub20deg <- preds_sa[["aligned_Slope"]] < 20
preds_notmire <- mask(preds_sa, mire, inverse = TRUE)
preds_notmiresub20 <- mask(preds_notmire, sub20deg, maskvalue = 0)
plot(preds_notmiresub20)
sum(complete.cases(values(preds_notmiresub20)))
set.seed(123)
osample <- eSample(envData = preds_notmiresub20, 
                   nExpect = 100, 
                   saveShape = FALSE,
                   nf = 2, 
                   lowerLim = 1e-2, upperLim = 1-1e-2)
ocoords <- st_as_sf(osample$GeoSamples)
ocoords <- st_set_crs(ocoords, 25833)
#st_write(ocoords, "ocoords.csv")

