library(tidyverse)
library(sf)
library(terra)
library(clhs)

N <- 160
cov <- rast("output/samplingCovariates.tif")

#Read data with coordinates and other attributes
df <- read_csv(file="output/preparedSamplingMatrix.csv")
vars <- c('RAD_K','RAD_TC','RAD_Th','RAD_U','slope')
summary(df)
cor(select(df, all_of(vars)))

dfcovs <- df[,c(vars, "rwalkcost")]
tictoc::tic()
clhsseeds <- 11 %>% 
  map(\(x) {
    set.seed(x)
    clhs(dfcovs, size = N, cost = "rwalkcost", iter = 1e6, simple = FALSE)
  })
tictoc::toc()
write_rds(clhsseeds, file="output/makeCLHS.rds")
# clhsseeds <- read_rds("output/makeCLHS.rds")

map_dbl(clhsseeds, \(x) x$obj[1e6]) %>% 
  plot()
ss <- clhsseeds[[which.min(map_dbl(clhsseeds, \(x) x$obj[1e6]))]] 

plot(ss); ss$obj[1e6]
s.dfcovs <- dfcovs[ss$index_samples,]
cor(select(s.dfcovs, all_of(vars)))
mySample <- df[ss$index_samples,]

sample <- st_as_sf(mySample, coords = c('x','y'), crs = crs(cov))
plot(cov['slope'], ext = ext(sample))
plot(st_geometry(sample), add=TRUE)
st_write(sample, "output/sampling.gpkg", layer = "cLHSseed15", append=FALSE)

# Monotone scaling of rwalkcost has no effect
# temp <- dfcovs
# temp <- mutate(temp, rwalkcost = rwalkcost^(2))
# set.seed(12)
# ss <- clhs(temp, size = N, cost = "rwalkcost")
# s.dfcovs <- dfcovs[ss,]
# cor(select(s.dfcovs, all_of(vars)))
# mySample <- df[ss,]

# clhs
# methods(clhs)
# clhs:::clhs.data.frame
