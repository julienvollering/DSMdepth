library(terra)
library(sf)
library(dplyr)
library(readr)
library(ggplot2)

rad.files <- list.files("data/", pattern = "RAD_miclev_.*10m.tif$",
                        full.names = TRUE)
rad10m <- rast(rad.files)
rad10m.mask <- st_read("data/Orskogfjellet-site.gpkg", layer="RAD_10m_mask")
rad10m <- mask(rad10m, rad10m.mask)
plot(rad10m)

slope1m <- rast("data/slope.tif")
slope1m

slope10m <- project(slope1m, rad10m, method = "average")
plot(slope10m)

covariates <- c(rad10m, slope10m)
names(covariates) <- gsub("_miclev", "", names(covariates))
names(covariates) <- gsub("_10m", "", names(covariates))

probe <- read_csv("data/depth/concatenatedRawProbe.csv")
probe <- st_as_sf(probe, coords = c('utm32E', 'utm32N'), crs = 25832)
probe <- select(probe, depth_cm) |> 
  st_transform(crs = crs(covariates))
df <- terra::extract(covariates, vect(probe), ID=FALSE) |> 
  bind_cols(probe)

df <- as_tibble(df) |> 
  select(-geometry)

pairs(df)
cor(df, use = "complete.obs")
corrplot::corrplot(cor(df, use = "complete.obs"), method = 'ellipse',diag=T,addCoef.col="black")

filter(df, !complete.cases(df))
df <- tidyr::drop_na(df)

dflong <- tidyr::pivot_longer(df, cols = 1:5, names_to = "covariate")
annot <- tribble(
  ~covariate, ~pearsonR,
  "RAD_K",    -0.23,
  "RAD_TC",   -0.21,
  "RAD_Th",   -0.29,
  "RAD_U",     0.00,
  "slope",    -0.28
)
ggplot(dflong) + 
  geom_point(aes(x = value, y = depth_cm)) + 
  facet_grid(cols = vars(covariate), scales = "free_x") +
  geom_text(
    data    = annot,
    mapping = aes(x = Inf, y = Inf, label = paste("Pearson r =",pearsonR)),
    hjust = 1.5, vjust = 10
  )
ggsave("output/exploratoryDataAnalysis.svg", height = 10, width = 5*10, unit="cm")
