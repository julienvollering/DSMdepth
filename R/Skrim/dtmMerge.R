library(terra) #v1.7-37
files <- list.files("data/Haram Skodje Ørskog Vestnes 2pkt 2015/data/dtm", pattern = ".tif$", full.names = TRUE)
tiles <- sprc(files)
summary(tiles)
merged <- merge(tiles)
merged
writeRaster(merged, "data/DTMmerged.tif")
