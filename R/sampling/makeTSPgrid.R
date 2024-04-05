library(tidyverse)
library(sf)
library(TSP)

# Illustration of TSP ####
set.seed(42)
x <- matrix(sample.int(10, 18, replace = TRUE), ncol = 2)
x
plot(x, pch=as.character(1:9), xlim=c(1,10), ylim=c(1,10))

dist(x)
dist(x, method="manhattan")

library(TSP)
?TSP
rtsp <- TSP(dist(x, method="manhattan"))
tour <- solve_TSP(rtsp)
tour

str(tour)

# TSP for making transects #### 
# copy attibute table from QGIS with selected features only
sfc <- read_delim(clipboard()) %>%
  pull(wkt_geom) %>% 
  st_as_sfc()
d <- read_delim(clipboard()) %>% 
  pull(for.cluster)
pts <- tibble(for.cluster = d, geom = sfc) %>% 
  st_as_sf()

rtsp <- TSP(dist(st_coordinates(pts), method="manhattan"), labels = pts$for.cluster)
rtsp <- insert_dummy(rtsp)
tour <- solve_TSP(rtsp)
tour
names(tour)
