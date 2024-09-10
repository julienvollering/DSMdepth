library(terra) #v1.7.78
files <- c("data/Skrim/eksport_926374_20240909/dtm1/data/dtm1_33_119_112.tif",
           "data/Skrim/eksport_926374_20240909/dtm1/data/dtm1_33_120_112.tif")
tiles <- sprc(files)
summary(tiles)
merged <- merge(tiles)
merged
writeRaster(merged, "data/Skrim/DTMmerged.tif", overwrite = TRUE)
