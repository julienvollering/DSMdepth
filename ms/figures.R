library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(sf)
library(rnaturalearth)

# Sites map ####

ne_download(scale = 50, type = 'populated_places', category = 'cultural', 
            destdir = "ms/figures/data", load = FALSE)
# ne_download(scale = 10, type = 'minor_islands', category = 'physical', 
#             destdir = here("data", "raw", "cartographic"), load = FALSE)
# 
cities <- ne_load(scale = 50, type = 'populated_places', category = 'cultural', 
                  destdir = "ms/figures/data", returnclass = 'sf')
norcities <- filter(cities, NAME %in% c("Oslo", "Trondheim", "Bergen"))
# norislands <- ne_load(scale = 10, type = 'minor_islands', category = 'physical', 
#                 destdir = here("data", "raw", "cartographic"), returnclass = 'sf') %>% 
#   st_crop(xmin = 4, xmax = 32, ymin = 58, ymax = 72) %>% 
#   st_transform(crs = 25833) # ETRS89 / UTM zone 33N

nor50 <- ne_countries(scale = 50, country = c("Norway"), returnclass = "sf") %>% 
  st_transform(crs = 25833)
orskog <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea")
skrim <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm")
both <- bind_rows(orskog,skrim)
labels <- both |> 
  st_centroid() |> 
  mutate(site = c("Ørskogfjellet", "Skrimfjella"))

CBdark2.3class <- c('#d95f02', '#7570b3', '#1b9e77')

g1 <- ggplot() +
  geom_sf(data = nor50, fill = "white") + 
  geom_sf(data = norcities, color = "grey50") +
  geom_sf_text(data = norcities, aes(label = NAME), color = "grey50",
               nudge_y = c(-1.5e4,1.5e4,2e4), 
               nudge_x = c(0, 4e4, 0), size = 2) +
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
    nudge_x = c(15e4, 5e4),
    segment.curvature = 0.5) +
  coord_sf(xlim = c(-7e4, 3.7e5), ylim = c(6.47e6, 7.05e6)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(g1, filename = 'sites-map.svg', path = "ms/figures",
       scale = 1, width = 74, height = 105, units = 'mm', dpi = 300)



# sessionInfo ####

sessioninfo::session_info()
