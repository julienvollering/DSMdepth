library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)

# Sites map ####

## Inset ####

# ne_download(scale = 50, type = 'populated_places', category = 'cultural', 
#             destdir = "ms/figures/data", load = FALSE)
cities <- ne_load(scale = 50, type = 'populated_places', category = 'cultural', 
                  destdir = "ms/figures/data", returnclass = 'sf')
norcities <- filter(cities, NAME %in% c("Oslo", "Trondheim", "Bergen"))

nor50 <- ne_countries(scale = 50, country = c("Norway"), returnclass = "sf") %>% 
  st_transform(crs = 25833)
orskog <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea") |> 
  st_transform(25833)
skrim <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_transform(25833)
both <- bind_rows(orskog, skrim)
labels <- both |> 
  st_centroid() |> 
  mutate(site = c("Ørskogfjellet", "Skrimfjella"))

g1 <- ggplot() +
  geom_sf(data = nor50, fill = "white") + 
  geom_sf(data = norcities, color = "grey50") +
  geom_sf_text(data = norcities, aes(label = NAME), color = "grey50",
               nudge_y = c(-2e4,1.5e4,2e4), 
               nudge_x = c(0, 5e4, 0), size = 2) +
  geom_sf(data = both,
          mapping = aes(color = 'red', fill = 'red')) +
  geom_text_repel(
    data = labels,
    aes(label = site, geometry = geom),
    stat = "sf_coordinates",
    min.segment.length = 0,
    arrow = arrow(length = unit(0.02, "npc")),
    force_pull = 0.1,
    nudge_y = c(-8e4, 12e4),
    nudge_x = c(13e4, 4e4),
    segment.curvature = 0.5) +
  coord_sf(xlim = c(-7e4, 3.7e5), ylim = c(6.47e6, 7.05e6), crs = 25833) +
  theme(legend.position = "none",
        axis.title.x = element_blank(), axis.title.y= element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt"))

## AR50 Ørskogfjellet ####

ar50aalesund <- st_read("ms/figures/data/1508_25832_ar50_gdb/1507_25832_ar50_gdb.gdb")
ar50vestnes <- st_read("ms/figures/data/1535_25832_ar50_gdb/1535_25832_ar50_gdb.gdb")
ar50orskog <- bind_rows(ar50aalesund,ar50vestnes) |> 
  mutate(artype = as.factor(artype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(geo)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)

# Creating hillshade
dtm <- rast("ms/figures/data/940695_dtm50/data/dtm50_6900_50m_33.tif")
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50orskog)
names(hill) <- "shades"
pal_greys <- hcl.colors(1000, "Grays")
pal_greys_index <- hill %>%
  mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)

# https://colorbrewer2.org/?type=qualitative&scheme=Paired&n=6
CBpaired.6class <- c('#1f78b4','#a6cee3','#33a02c','#b2df8a','#e31a1c','#fb9a99')

g2 <- ggplot(ar50orskog) +
  geom_spatraster(data = hill, fill = pal_greys[pal_greys_index], maxcell = Inf, 
                  alpha = 1) +
  geom_sf(data = ar50orskog, 
          mapping = aes(fill = artype), alpha = 0.6, color = NA) +
  geom_sf(data = orskog, fill = NA, color = 'black', 
          linewidth = 1, linetype = 1) + 
  coord_sf(xlim = st_bbox(orskog)[c(1,3)], ylim = st_bbox(orskog)[c(2,4)],
           crs = 25833) +
  scale_fill_discrete(labels = c("Built-up",
                                 "Agricultural",
                                 "Forest",
                                 "Open upland",
                                 "Peatland",
                                 "Freshwater"), type = CBpaired.6class) +
  guides(fill = guide_legend(title = NULL)) +
  annotation_scale(location = "tl") + 
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.margin = margin(1, 1, 1, 1, "pt")
  )

## AR50 Skrim ####

ar50kongsberg <- st_read("ms/figures/data/3303_25833_ar50_gml/3006_25833_ar50_gml.gml",
                         "ArealressursFlate") |> 
  st_set_crs(25833)
ar50skrim <- ar50kongsberg |> 
  mutate(artype = as.factor(arealtype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(område)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)

# Creating hillshade
dtm <- c("ms/figures/data/940334_dtm50/data/dtm50_6601_50m_33.tif",
         "ms/figures/data/940334_dtm50/data/dtm50_6602_50m_33.tif") |> 
  sprc() |> 
  merge()
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50skrim)
names(hill) <- "shades"
pal_greys <- hcl.colors(1000, "Grays")
pal_greys_index <- hill %>%
  mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)

g3 <- ggplot(ar50skrim) +
  geom_spatraster(data = hill, fill = pal_greys[pal_greys_index], maxcell = Inf, 
                  alpha = 1) +
  geom_sf(data = ar50skrim, 
          mapping = aes(fill = artype), alpha = 0.6, color = NA) +
  geom_sf(data = skrim, fill = NA, color = 'black', 
          linewidth = 1, linetype = 1) + 
  coord_sf(xlim = st_bbox(skrim)[c(1,3)], ylim = st_bbox(skrim)[c(2,4)],
           crs = 25833) +
  scale_fill_discrete(type = CBpaired.6class) +
  guides(fill = "none") +
  annotation_scale(location = "bl") +
  theme(plot.margin = margin(1, 1, 1, 1, "pt"))

## Combined ####

# Patchwork nesting
# g2 /(g3 | g1) #too much whitespace (horizontally between but especially and unnecessarily on sides)
# g2 | (g3 / g1) #too much whitespace (vertically, but not on top/bottom)
# seems to be an issue with aligning geom_sf (fixed-aspect) plots

#saving to device
layout <- c( #use 10x10 layout to adjust design
  patchwork::area(1,1,7,10),
  patchwork::area(8,1,10,7),
  patchwork::area(8,8,10,10)
)
plot(layout)
g2 + g3 + g1 + plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'sites-patchwork.pdf', path = "ms/figures",
       width =210-30, height = 240-40, units = 'mm') #copernicus.cls page 210x240

## With tmap instead of ggplot2 ####
# library(tmap)
# library(tmaptools)
# 
# tm_shape(ar50, bbox=orskog) +
#   tm_polygons("artype") +
#   tm_shape(orskog) +
#   tm_borders(col = 'red')
# 
# # With background tiles from OSM
# osm_orskog <- read_osm(orskog, ext=1.1,
#                        type = "osm")
# tm_shape(osm_orskog) +
#   tm_rgb() +
# tm_shape(orskog) +
#   tm_borders()
# 
# map <- OpenStreetMap::openmap(c(62,7),
#                              c(61.5,7.5),
#                              type='osm')
# map <- OpenStreetMap::openmap(c(62,7),
#                               c(61.5,7.5),
#                               type='https://opencache.statkart.no/gatekeeper/gk/gk.open_gmaps?layers=topo4&zoom={z}&x={x}&y={y}')
# plot(map)

# sessionInfo ####

sessioninfo::session_info()
