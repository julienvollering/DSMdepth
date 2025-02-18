library(tidyverse)
library(sf)

peatDataOrder <- c("recordDate", "surveyor", "institution", "location", 
                   "country", "latitude", "longitude", "coordinateUncertainty",
                   "geodeticDatum", "gridReference", "easting", "northing",
                   "locationRemarks", "surfacePeatDepth", "peatDepthResolution",
                   "peatPresence", "buriedPeatStartDepth", "buriedPeatEndDepth",
                   "surfaceElevation", "probeType", "probeReachedBottom", 
                   "peatDepthRemarks")

# Ørskogfjellet ####

## Probes ####

probe <- read_csv("data/depth/concatenatedCleanProbe.csv")
standard <- st_as_sf(probe, coords = c('X', 'Y', 'utm32Z'), crs = 25832) |> 
  st_transform(4326) |> 
  select(HRMS, depth_cm, comment)
standard <- standard |> 
  mutate(probeReachedBottom = if_else(
    grepl("plus|\\+", comment, ignore.case = TRUE), 
    "No", "Yes"))

filter(standard, !is.na(comment)) |> 
  pull(comment) |> 
  unique()
standard <- standard |> 
  mutate(
    peatDepthRemarks = 
      case_when(
        comment == "gravemasser veg" ~ "fill material road",                               
        comment == "fyllmasse" ~ "fill material",                                     
        comment == "edge of infill" ~ comment,                                
        comment == "stein" ~ "stone",                                         
        comment == "tue" ~ "peat hummock",                                           
        comment == "+, naadde ikke bunn, tue" ~ "peat hummock",                      
        comment == "fyllmasse i  veggroft" ~ "fill material in road ditch",                         
        comment == "nederste 30 cm er sikker sand, muligens hoyere" ~ "bottom 30 cm may be sand, possibly higher",  
        comment == "ikke tydelig bunn" ~ "no clear bottom",                             
        comment == "uklar grense, fin sand i randen paa torvspydet" ~ "unclear boundary, fine sand in the edge of the peat probe",
           ),
    locationRemarks =
      case_when(
        comment == "float  tett skog" ~ "dense forest",                              
        comment == "next to ditch, float" ~ "next to ditch",                     
        comment == "slope towards ditch, float" ~ "slope towards ditch",         
        comment == "in ditch, float" ~ "in ditch",
        comment == "in wide running ditch" ~ comment,                
        comment == "in ditch" ~ comment,
        comment == "outside ditch" ~ comment,                                 
        comment == "road edge disturbed" ~ comment,                           
        comment == "road" ~ comment,                                          
        comment == "next to road" ~ comment,                                  
        comment == "in stream" ~ comment,                                     
        comment == "next to stream" ~ comment,                               
        comment == "between road and stream" ~ comment,                       
        comment == "close to stream" ~ comment,                               
        comment == "ditch" ~ comment,                                         
        comment == "edge ditch" ~ comment,
        comment == "float skog" ~ "forest"
      ))
filter(standard, grepl("fill", peatDepthRemarks, ignore.case = TRUE) | 
         grepl("fill", locationRemarks, ignore.case = TRUE))

standard <- standard |> 
  mutate(coordinateUncertainty = HRMS,
         surfacePeatDepth = depth_cm) |> 
  select(-HRMS, -depth_cm, -comment) 

standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y, surfaceElevation = Z,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "2023-08",
         surveyor = "Karl-Kristian Muggerud | Sigurd Daniel Nerhus | Knut Rydgren | Julien Vollering",
         institution = "Western Norway University of Applied Sciences",
         location = "Møre og Romsdal County",
         country = "Norway",
         peatDepthResolution = "3 cm",
         probeType = "peat probe")

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/OrskogfjelletProbe.csv", append = FALSE)

## GPR ####

gpr <- read_csv("data/depth/all.csv") |> 
  filter(source == "gpr") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

standard <- gpr |> 
  st_transform(4326) |> 
  select(surfacePeatDepth = depth_cm)
  
standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "2023-08",
         surveyor = "Karl-Kristian Muggerud | Sigurd Daniel Nerhus | Knut Rydgren | Julien Vollering",
         institution = "Western Norway University of Applied Sciences",
         location = "Møre og Romsdal County",
         country = "Norway",
         coordinateUncertainty = "30", #https://dwc.tdwg.org/list/#dwc_coordinateUncertaintyInMeters
         peatDepthResolution = "42 cm", # MAE in wave velocity calibration
         probeType = "ground-penetrating radar",
         probeReachedBottom = "Yes",
         peatDepthRemarks = "peatDepthResolution is set to the mean absolute error obtained in a wave velocity calibration (which also includes error from lateral separation) and we believe measurement accuracy to be better in most measurements."
         )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/OrskogfjelletGPR.csv", append = FALSE)

## Wisen2021 ####

wisen <- read_csv("data/depth/all.csv") |> 
  filter(source == "wisen2021") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

standard <- wisen |> 
  st_transform(4326) |> 
  select(surfacePeatDepth = depth_cm)

standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "2020/2021",
         surveyor = "Roger Wisén | Fredrik Olsen | IMPAKT Geofysik AB",
         institution = "The Norwegian Public Roads Administration",
         location = "Møre og Romsdal County",
         country = "Norway",
         probeType = "ground-penetrating radar",
         probeReachedBottom = "Yes"
  )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/OrskogfjelletNorwegianPublicRoadsAdministration.csv", append = FALSE)

## Myrarkivet ####

myrarkivet <- read_csv("data/depth/all.csv") |> 
  filter(source == "myrarkivet") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

standard <- myrarkivet |> 
  st_transform(4326) |> 
  select(surfacePeatDepth = depth_cm)

standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "1983-07-11",
         surveyor = "The Norwegian Soil and Mire Company",
         institution = "Norwegian Institute of Bioeconomy Research",
         location = "Møre og Romsdal County",
         country = "Norway",
         peatDepthResolution = "10 cm",
         probeReachedBottom = "Yes"
  )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/OrskogfjelletNorwegianInstituteBioeconomyResearch.csv", append = FALSE)

# Skrimfjella ####

## Probes ####

probe <- read_csv("data/Skrim/depth_probe.csv")
standard <- st_as_sf(probe, coords = c('X', 'Y'), crs = 4326) |> 
  select(depth_cm, note)

filter(standard, !is.na(note)) |> 
  pull(note) |> 
  unique()

standard <- standard |> 
  mutate(probeReachedBottom = if_else(
    grepl("plus|\\+", note, ignore.case = TRUE) & depth_cm == 570, 
    "No", "Yes"))

standard <- standard |> 
  mutate(
    locationRemarks =
      case_when(
        note == "River" ~ note,                              
        note == "Forest" ~ note,                              
        note == "Hillside" ~ note                              
      ))

standard <- standard |> 
  mutate(coordinateUncertainty = 5, # 3 m GPS + 2 m triangle
         surfacePeatDepth = depth_cm) |> 
  select(-depth_cm, -note) 

standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "2020-08",
         surveyor = "Karl-Kristian Muggerud | Mikko Sparf | Mette Kusk Gillespie | Knut Rydgren",
         institution = "Western Norway University of Applied Sciences",
         location = "Buskerud County",
         country = "Norway",
         peatDepthResolution = "3 cm",
         probeType = "peat probe"
         )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/SkrimfjellaProbe.csv", append = FALSE)

## GPR ####

gpr <- read_csv("data/Skrim/depth_all.csv") |> 
  filter(source == "gpr") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25833)

standard <- gpr |> 
  st_transform(4326) |> 
  select(surfacePeatDepth = depth_cm)

standard <- standard |> 
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused")

standard <- standard |>
  mutate(recordDate = "2020-08",
         surveyor = "Karl-Kristian Muggerud | Mikko Sparf | Mette Kusk Gillespie | Knut Rydgren",
         institution = "Western Norway University of Applied Sciences",
         location = "Buskerud County",
         country = "Norway",
         coordinateUncertainty = "30", #https://dwc.tdwg.org/list/#dwc_coordinateUncertaintyInMeters
         peatDepthResolution = "29 cm", # MAE in wave velocity calibration
         probeType = "ground-penetrating radar",
         probeReachedBottom = "Yes",
         peatDepthRemarks = "peatDepthResolution is set to the mean absolute error obtained in a wave velocity calibration (which also includes error from lateral separation) and we believe measurement accuracy to be better in most measurements."
  )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/SkrimfjellaGPR.csv", append = FALSE)

## Occurrence ####

occurrence <- st_read("data/Skrim/Skrim-site.gpkg", "occurrence")

standard <- occurrence |> 
  st_transform(4326) |> 
  mutate(peatPresence = if_else(peat == 1, "Yes", "No")) |> 
  select(peatPresence)

standard <- standard |>
  bind_cols(st_coordinates(standard)) |> 
  st_drop_geometry() |> 
  mutate(longitude = X, latitude = Y,
         geodeticDatum = "https://epsg.io/4326", .keep = "unused") |> 
  as_tibble()

standard <- standard |>
  mutate(recordDate = "2020-08",
         surveyor = "Karl-Kristian Muggerud | Mikko Sparf",
         institution = "Western Norway University of Applied Sciences",
         location = "Buskerud County",
         country = "Norway",
         coordinateUncertainty = "5", 
         probeType = "other",
         probeReachedBottom = "No",
         peatDepthRemarks = "peatPresence was determined by inspecting the top 20 cm of soil with a trowel."
  )

standard <- standard |> 
  select(all_of(intersect(peatDataOrder, names(standard))))

write_csv(standard, "data/depth/archive/SkrimfjellaOccurrence.csv", append = FALSE)

# Metadata ####

metadataOrder <- c(
  "title", "organizationName", "abstract", "additionalInfo", "license", "rights",
  "bibliographicCitation", "geographicDescription", "westBoundingCoordinate", 
  "eastBoundingCoordinate", "northBoundingCoordinate", "southBoundingCoordinate",
  "temporalCoverageStartDate", "temporalCoverageEndDate", "purpose", "methods",
  "qualityControl", "contactGivenName", "contactSurname", "contactElectronicMailAddress",
  "contactPhone"
  )

metadata <- list.files("data/depth/archive", pattern = ".csv") |> 
  map(~ read_csv(paste0("data/depth/archive/", .))) |> 
  set_names(list.files("data/depth/archive", pattern = ".csv")) |> 
  map(~ mutate(., across(recordDate, as.character))) |>
  list_rbind(names_to = "file")

standard <- metadata |> 
  group_by(file) |> 
  summarize(
    westBoundingCoordinate = min(longitude),
    eastBoundingCoordinate = max(longitude),
    northBoundingCoordinate = max(latitude),
    southBoundingCoordinate = min(latitude)
  )

joinMetadata <- readxl::read_xlsx("data/depth/archive/joinMetadata.xlsx")

standard <- standard |> 
  left_join(joinMetadata, by = c("file" = "file"))

standard <- standard |> 
  select(all_of(c("file", intersect(metadataOrder, names(standard)))))

write_csv(standard, "data/depth/archive/fileMetadata.csv", append = FALSE)

siteCoverage <- metadata |> 
  group_by(location) |> 
  summarize(
    westBoundingCoordinate = min(longitude),
    eastBoundingCoordinate = max(longitude),
    northBoundingCoordinate = max(latitude),
    southBoundingCoordinate = min(latitude)
  )

