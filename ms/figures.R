# Sites map ####

library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)

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
  mutate(site = c("\u00D8rskogfjellet", "Skrimfjella"))

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
# ggsave(filename = 'sites-patchwork.svg', path = "ms/figures",
#        width =210-30, height = 240-40, units = 'mm')

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

# Model metrics ####

library(tidyverse)

orskog <- read_csv("output/modelmetrics.csv") 
skrim <- read_csv("output/Skrim/modelmetrics.csv")
plotting <- bind_rows(orskog = orskog, skrim = skrim, .id = 'site') %>% 
  filter(.metric != 'mae') %>% 
  mutate(
    site = fct_relevel(site, "orskog"),
    .metric = fct_relevel(.metric, "ccc", "rsq", "rmse"),
    metric = case_when(
      .metric == "ccc" ~ "Concordance \ncorrelation coefficient",
      .metric == "rsq" ~ "R-squared",
      .metric == "rmse" ~ "RMSE"),
    model = case_when(
      model == "DMK" ~ "DMK class (2)",
      model == "Terrain" ~ "terrain (21)",
      model == "TerrainDMK" ~ "terrain + DMK class (23)",
      model == "RadiometricTerrain" ~ "terrain + radiometric (25)",
      model == "RadiometricTerrainDMK" ~ "all predictors (27)"),
    model = fct_relevel(model,
                        "DMK class (2)",
                        "terrain (21)",
                        "terrain + DMK class (23)",
                        "terrain + radiometric (25)",
                        "all predictors (27)"))
str(plotting)

g1 <- ggplot(plotting) +
  geom_linerange(aes(y = model, 
                     x = mean, 
                     xmin = mean - std_err, 
                     xmax = mean + std_err, 
                     color = site),
                 position = position_dodge2(width= 0.3, reverse = TRUE)) +
  geom_point(aes(y = model, x = mean, color = site),
             position = position_dodge2(width= 0.3, reverse = TRUE)) +
  geom_text(aes(y = model, x = mean, label = signif(mean, 2), group = site),
            position = position_dodge2(width= 1, reverse = TRUE),
            color = "grey50", size = 2.5) +
  facet_wrap(~metric, nrow = 1, scales = "free_x", axes ="margins") +
  scale_color_discrete(labels = c(orskog = "\u00D8rskogfjellet", 
                                  skrim = "Skrimfjella")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(-0.15, 1.08),
        legend.key.spacing.y = unit(-2, "mm"),)
ggsave(filename = 'modelmetrics.pdf', path = "ms/figures",
       width =210-30, height = (240-40)/2.5, units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'modelmetrics.svg', path = "ms/figures",
#        width =210-30, height = (240-40)/2.5, units = 'mm')

# Variable importance ####

library(tidyverse)
library(patchwork)

orskog <- read_csv("output/variable_importance.csv") %>% 
  filter(type != "perm.ranger")
skrim <- read_csv("output/Skrim/variable_importance.csv") %>% 
  filter(type != "perm.ranger")

orskog.rank <- orskog %>% 
  group_by(Variable) %>% 
  summarize(Imp.median = median(Importance))
skrim.rank <- skrim %>% 
  group_by(Variable) %>% 
  summarize(Imp.median = median(Importance))

plot.orskog <- left_join(orskog, orskog.rank, by = join_by(Variable))
labs.orskog <- plot.orskog |> 
  distinct(Variable, removed, Imp.median) |>
  replace_na(list(removed = "")) |>
  arrange(Imp.median) |> 
  pull(removed) |> 
  stringr::str_wrap(width = 40)
g1 <- plot.orskog %>% 
  mutate(Variable = case_when(
    Variable == "dmkdepth_grunn.myr" ~ "DMK_shallow",
    Variable == "dmkdepth_unknown" ~ "DMK_unknown",
    TRUE ~ Variable)) %>%
  ggplot() +
  geom_col(aes(x = fct_reorder(Variable, Imp.median), 
               y = Importance,
               group = type,
               fill = type),
           position = position_dodge2()) +
  scale_x_discrete(sec.axis = dup_axis(labels = labs.orskog)) +
  coord_flip() +
  scale_fill_discrete(labels = c(firm = "FIRM", 
                                  perm.vip = "permutation", 
                                  shap = "Shapley")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y.right = element_text(size = 8, color = "grey50"),
        legend.position = "inside",
        legend.position.inside = c(1.2, 0.85)) + 
  guides(fill = guide_legend(reverse = TRUE))

plot.skrim <- left_join(skrim, skrim.rank, by = join_by(Variable))
labs.skrim <- plot.skrim |> 
  distinct(Variable, removed, Imp.median) |>
  replace_na(list(removed = "")) |>
  arrange(Imp.median) |> 
  pull(removed) |> 
  stringr::str_wrap(width = 40)
g2 <- plot.skrim %>% 
  mutate(Variable = case_when(
    Variable == "dmkdepth_grunn" ~ "DMK_shallow",
    Variable == "dmkdepth_unknown" ~ "DMK_unknown",
    TRUE ~ Variable)) %>%
  ggplot() +
  geom_col(aes(x = fct_reorder(Variable, Imp.median), 
               y = Importance,
               group = type,
               fill = type),
           position = position_dodge2()) +
  scale_x_discrete(sec.axis = dup_axis(labels = labs.skrim)) +
  coord_flip() +
  guides(fill = "none") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y.right = element_text(size = 8, color = "grey50"))

## Combined ####

g1 + g2 + plot_layout(ncol = 1, axes = 'collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'variable_importance.pdf', path = "ms/figures",
       width =(210-30), height = (240-40)/1.25, units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'variable_importance.svg', path = "ms/figures",
#        width =(210-30)/1.5, height = (240-40)/1.25, units = 'mm')

# Partial dependence plots ####

library(tidyverse)
library(patchwork)

orskog <- read_csv("output/pdpice.csv")
orskog.vi <- read_csv("output/variable_importance.csv") %>% 
  filter(type != "perm.ranger")
orskog.rank <- orskog.vi %>% 
  group_by(Variable) %>% 
  summarize(Imp.median = median(Importance)) %>% 
  arrange(desc(Imp.median)) %>% 
  slice_head(n = 6)
orskog <- left_join(orskog, orskog.rank, by = join_by(feature == Variable)) %>%
  drop_na(Imp.median) %>%
  mutate(
    feature = case_when(
      feature == "dmkdepth_grunn.myr" ~ "DMK_shallow",
      feature == "dmkdepth_unknown" ~ "DMK_unknown",
      TRUE ~ feature),
    feature = fct_reorder(feature, -Imp.median))

g1 <- ggplot(data = orskog, aes(x = .borders, y = .value, group = .id)) + 
  geom_line(data = filter(orskog, .type == 'ice'), alpha = 0.01) +
  geom_line(data = filter(orskog, .type == 'pdp'), color = "red") +
  geom_rug(data = filter(orskog, .type == "observed"), sides = "b") +
  coord_cartesian(ylim = c(0, 400)) +
  facet_wrap(~feature, scales = "free_x") + 
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank())

skrim <- read_csv("output/Skrim/pdpice.csv")
skrim.vi <- read_csv("output/Skrim/variable_importance.csv") %>% 
  filter(type != "perm.ranger")
skrim.rank <- skrim.vi %>% 
  group_by(Variable) %>% 
  summarize(Imp.median = median(Importance)) %>% 
  arrange(desc(Imp.median)) %>% 
  slice_head(n = 6)
skrim <- left_join(skrim, skrim.rank, by = join_by(feature == Variable)) %>%
  drop_na(Imp.median) %>%
  mutate(
    feature = case_when(
      feature == "dmkdepth_grunn.myr" ~ "DMK_shallow",
      feature == "dmkdepth_unknown" ~ "DMK_unknown",
      TRUE ~ feature),
    feature = fct_reorder(feature, -Imp.median))

g2 <- ggplot(data = skrim, aes(x = .borders, y = .value, group = .id)) + 
  geom_line(data = filter(skrim, .type == 'ice'), alpha = 0.01) +
  geom_line(data = filter(skrim, .type == 'pdp'), color = "red") +
  geom_rug(data = filter(skrim, .type == "observed"), sides = "b") +
  coord_cartesian(ylim = c(0, 250)) +
  facet_wrap(~feature, scales = "free_x") + 
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank())

g1 + g2 + plot_layout(ncol = 1, guides = 'collect', axes = 'collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'partial_dependence.pdf', path = "ms/figures",
       width =(210-30), height = (240-40), units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'partial_dependence.svg', path = "ms/figures",
#        width =(210-30), height = (240-40), units = 'mm')
  
# sessionInfo ####

sessioninfo::session_info()
