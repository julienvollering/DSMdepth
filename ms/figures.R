# Sites map 2-panel ####

library(rnaturalearth)
library(sf)
library(terra)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggspatial)
library(ggnewscale)
library(scales)

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

pNorway <- ggplot() +
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
    arrow = arrow(length = unit(0.05, "npc")),
    direction = "y",
    force_pull = 0.1,
    nudge_y = c(-10e4, 12e4),
    nudge_x = c(3e4, 4e4),
    segment.curvature = 0.2) +
  coord_sf(xlim = c(-6.1e4, 3.7e5), ylim = c(6.47e6, 7.05e6), crs = 25833) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y= element_blank(),
    axis.ticks = element_blank(), axis.text = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

## Ørskogfjellet ####

orskog <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea") |> 
  st_transform(25833)
bboxOrskog <- st_bbox(st_buffer(orskog, 500))

ar50aalesund <- st_read("ms/figures/data/1508_25832_ar50_gdb/1507_25832_ar50_gdb.gdb")
ar50vestnes <- st_read("ms/figures/data/1535_25832_ar50_gdb/1535_25832_ar50_gdb.gdb")
ar50orskog <- bind_rows(ar50aalesund,ar50vestnes) |> 
  mutate(artype = as.factor(artype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(geo)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)
ar50orskogmask <- ar50orskog |> 
  group_by(TRUE) %>%
  summarize()

# Creating hillshade
dtm <- rast("ms/figures/data/940695_dtm50/data/dtm50_6900_50m_33.tif") |> 
  crop(bboxOrskog)
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50orskogmask)
hilldf_singleOrskog <- as.data.frame(hill, xy = TRUE)

# Creating contours
elev_range <- dtm |> 
  values() |> 
  range(na.rm = TRUE)
contourlinesOrskog <- dtm |>
  stars::st_as_stars() |> 
  stars::st_contour(contour_lines = TRUE,
                    breaks = seq(from = floor(elev_range[1]), 
                                 to = ceiling(elev_range[2]), 
                                 by = 100))

# https://colorbrewer2.org/?type=qualitative&scheme=Paired&n=6
CBpaired.6class <- c('#e31a1c','#fb9a99','#33a02c','#b2df8a','#a6cee3','#1f78b4')

mapOrskog <- ggplot(ar50orskog) +
  geom_tile(data = hilldf_singleOrskog, 
            aes(x, y, fill = hillshade), show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_sf(data = ar50orskog,
          mapping = aes(fill = artype), alpha = 0.75, color = NA) +
  scale_fill_discrete(labels = c("Built-up",
                                 "Agricultural",
                                 "Forest",
                                 "Open upland",
                                 "Peatland",
                                 "Freshwater"), type = CBpaired.6class) +
  geom_sf(data = contourlinesOrskog, linewidth = 0.2, color = "#fdbf6f",
          show.legend = FALSE) +
  geom_sf(data = orskog,
          mapping = aes(color = 'Study area   '), fill = NA,
          linewidth = 0.5, linetype = 1) +
  scale_color_manual(values = 'black') +
  coord_sf(
    xlim = bboxOrskog[c(1,3)], 
    ylim = bboxOrskog[c(2,4)],
    crs = 25833, expand = FALSE, label_axes = "EN--") +
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL)) +
  annotation_scale(location = "tl", width_hint = 0.15) +
  annotation_north_arrow(
    location = "tl", 
    pad_x = unit(0, "cm"),
    pad_y = unit(1, "cm"),
    which_north = "true", 
    style = north_arrow_minimal()) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_blank(), axis.title.y= element_blank()
  )

pOrskog <- mapOrskog + labs(tag = "(b)")

## Skrim ####

skrim <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_transform(25833)
bboxSkrim <- st_bbox(st_buffer(skrim, 500))

ar50kongsberg <- st_read("ms/figures/data/3303_25833_ar50_gml/3006_25833_ar50_gml.gml",
                         "ArealressursFlate") |> 
  st_set_crs(25833)
ar50skrim <- ar50kongsberg |> 
  mutate(artype = as.factor(arealtype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(område)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)
ar50skrimmask <- ar50skrim |> 
  group_by(TRUE) %>%
  summarize()

# Creating hillshade
dtm <- c("ms/figures/data/940334_dtm50/data/dtm50_6601_50m_33.tif",
         "ms/figures/data/940334_dtm50/data/dtm50_6602_50m_33.tif") |> 
  sprc() |> 
  merge() |> 
  crop(bboxSkrim)
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50skrimmask)
hilldf_singleSkrim <- as.data.frame(hill, xy = TRUE)

# Creating contours
elev_range <- dtm |> 
  values() |> 
  range(na.rm = TRUE)
contourlinesSkrim <- dtm |>
  stars::st_as_stars() |> 
  stars::st_contour(contour_lines = TRUE,
                    breaks = seq(from = floor(elev_range[1]), 
                                 to = ceiling(elev_range[2]), 
                                 by = 100))

mapSkrim <- ggplot(ar50skrim) +
  geom_tile(data = hilldf_singleSkrim, 
            aes(x, y, fill = hillshade), show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_sf(data = ar50skrim,
          mapping = aes(fill = artype), alpha = 0.75, color = NA) +
  scale_fill_discrete(type = CBpaired.6class) +
  geom_sf(data = contourlinesSkrim, linewidth = 0.2, color = "#fdbf6f",
          show.legend = FALSE) +
  geom_sf(data = skrim, fill = NA, color = 'black', 
          linewidth = 0.5, linetype = 1) + 
  scale_color_manual(values = 'black') +
  coord_sf(
    xlim = bboxSkrim[c(1,3)], 
    ylim = bboxSkrim[c(2,4)],
    crs = 25833, expand = FALSE, label_axes = "EN--") +
  guides(fill = "none", color = "none") +
  annotation_scale(location = "tr") +
  annotation_north_arrow(
    location = "tr", 
    pad_x = unit(0, "cm"),
    pad_y = unit(1, "cm"),
    which_north = "true", 
    style = north_arrow_minimal()) +
  theme(
    axis.title.x = element_blank(), axis.title.y= element_blank()
  )

pSkrim <- mapSkrim + labs(tag = "(a)") +
  inset_element(pNorway, 
                left = 0, bottom = 0, right = 0.47, top = 0.78, 
                align_to = "panel", clip = TRUE, ignore_tag = TRUE)

## Combined ####

pAll <- pSkrim + pOrskog + plot_layout(ncol = 1, nrow = 2, heights = c(2, 3))
ggsave(pAll, filename = 'map-sites.pdf', path = "ms/figures",
       width =210-30, height = 297-30-40, units = 'mm') #A4 page 210x297

# Distribution map 2-panel ####

library(sf)
library(terra)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggspatial)
library(ggnewscale)
library(scales)

probs <- c(0.25, 0.5, 0.75)
binwidth <- 20

modelOrskog <- sf::st_read("output/modeling.gpkg", layer="dataframe") |> 
  filter(!(ar5cover %in% c(11,12)))
modelSkrim <- sf::st_read("output/Skrim/modeling.gpkg", layer="dataframe")|> 
  filter(!(ar5cover %in% c(11,12)))
color_limits <- range(c(modelOrskog$depth_cm, modelSkrim$depth_cm))

## Ørskogfjellet ####

orskog <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea") |> 
  st_transform(25833)
bboxOrskog <- st_bbox(st_buffer(orskog, 500))

ar50aalesund <- st_read("ms/figures/data/1508_25832_ar50_gdb/1507_25832_ar50_gdb.gdb")
ar50vestnes <- st_read("ms/figures/data/1535_25832_ar50_gdb/1535_25832_ar50_gdb.gdb")
ar50orskog <- bind_rows(ar50aalesund,ar50vestnes) |> 
  mutate(artype = as.factor(artype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(geo)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)
ar50orskogmask <- ar50orskog |> 
  group_by(TRUE) %>%
  summarize()

# Creating hillshade
dtm <- rast("ms/figures/data/940695_dtm50/data/dtm50_6900_50m_33.tif") |> 
  crop(bboxOrskog)
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50orskogmask)
hilldf_singleOrskog <- as.data.frame(hill, xy = TRUE)

mapOrskog <- ggplot(orskog) +
  geom_tile(data = hilldf_singleOrskog, 
            aes(x, y, fill = hillshade), show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  geom_sf(data = modelOrskog, 
          mapping = aes(color = depth_cm), alpha = 1) +
  scale_color_viridis_c(limits = color_limits, direction = -1, guide = "none") +
  geom_sf(data = orskog, color = 'black', fill = NA, 
          linewidth = 0.5, linetype = 1,
          show.legend = FALSE) +
  coord_sf(
    xlim = bboxOrskog[c(1,3)], 
    ylim = bboxOrskog[c(2,4)],
    crs = 25833, expand = FALSE, label_axes = "EN--") +
  annotation_scale(location = "tl", width_hint = 0.15) +
  annotation_north_arrow(
    location = "tl", 
    pad_x = unit(0, "cm"),
    pad_y = unit(1, "cm"),
    which_north = "true", 
    style = north_arrow_minimal()) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_blank(), axis.title.y= element_blank()
  )

### Depth distribution ####
mean_depth <- mean(modelOrskog$depth_cm)
quantiles_depth <- quantile(modelOrskog$depth_cm, c(0.25, 0.5, 0.75))

p_base <- ggplot(modelOrskog, aes(x = depth_cm)) +
  geom_histogram(binwidth = binwidth)
hist_data <- layer_data(p_base)
x_values <- c(mean_depth, quantiles_depth)
bin_indices <- findInterval(x_values, hist_data$xmin, rightmost.closed = TRUE)
y_values <- hist_data$count[bin_indices]

label_data <- data.frame(
  x = x_values,
  y = y_values,
  label = c(paste("mean:", round(mean_depth), "cm"),
            paste("Q1:", round(quantiles_depth[1]), "cm"),
            paste("median:", round(quantiles_depth[2]), "cm"),
            paste("Q3:", round(quantiles_depth[3]), "cm")
))

densityOrskog <- ggplot(modelOrskog, aes(x = depth_cm)) +
  geom_histogram(aes(fill = after_stat(x)), 
                 binwidth = binwidth, color = "black", size = 0.1) +
  scale_fill_viridis_c(limits = color_limits, direction = -1, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(modelOrskog$depth_cm), by = 200),
                     labels = scales::label_number(suffix = " cm")) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text_repel(data = label_data, aes(x = x, y = y, label = label),
                  direction = "y",
                  min.segment.length = 0,
                  nudge_x = max(hist_data$x) * 0.2,
                  nudge_y = max(hist_data$count) * 0.1,
                  size = 3) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.line.x = element_line()
  )

pOrskog <- mapOrskog + labs(tag = "(b)") +
  inset_element(densityOrskog, left = 0.47, bottom = 0.02, right = 0.89, top = 0.5, 
                align_to = "full", clip = TRUE, ignore_tag = TRUE)

## Skrim ####

skrim <- st_read("data/Skrim/Skrim-site.gpkg", "fieldsite_outline_utm") |> 
  st_transform(25833)
bboxSkrim <- st_bbox(st_buffer(skrim, 500))

ar50kongsberg <- st_read("ms/figures/data/3303_25833_ar50_gml/3006_25833_ar50_gml.gml",
                         "ArealressursFlate") |> 
  st_set_crs(25833)
ar50skrim <- ar50kongsberg |> 
  mutate(artype = as.factor(arealtype)) |> 
  group_by(artype) %>%
  summarize(geometry = st_union(område)) |> 
  filter(!(artype %in% c(70,82,99))) |> 
  st_transform(crs = 25833)
ar50skrimmask <- ar50skrim |> 
  group_by(TRUE) %>%
  summarize()

# Creating hillshade
dtm <- c("ms/figures/data/940334_dtm50/data/dtm50_6601_50m_33.tif",
         "ms/figures/data/940334_dtm50/data/dtm50_6602_50m_33.tif") |> 
  sprc() |> 
  merge() |> 
  crop(bboxSkrim)
slope <- terrain(dtm, "slope", unit = "radians")
aspect <- terrain(dtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 225)
hill <- mask(hill, ar50skrimmask)
hilldf_singleSkrim <- as.data.frame(hill, xy = TRUE)

mapSkrim <- ggplot(skrim) +
  geom_tile(data = hilldf_singleSkrim, 
            aes(x, y, fill = hillshade), show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  geom_sf(data = modelSkrim, 
          mapping = aes(color = depth_cm), alpha = 1) +
  scale_color_viridis_c(limits = color_limits, direction = -1, guide = "none") +
  geom_sf(data = skrim, color = 'black', fill = NA, 
          linewidth = 0.5, linetype = 1,
          show.legend = FALSE) +
  coord_sf(
    xlim = bboxSkrim[c(1,3)], 
    ylim = bboxSkrim[c(2,4)],
    crs = 25833, expand = FALSE, label_axes = "EN--") +
  annotation_scale(location = "tr") +
  annotation_north_arrow(
    location = "tr", 
    pad_x = unit(0, "cm"),
    pad_y = unit(1, "cm"),
    which_north = "true", 
    style = north_arrow_minimal()) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_blank(), axis.title.y= element_blank()
  )

### Depth distribution ####
mean_depth <- mean(modelSkrim$depth_cm)
quantiles_depth <- quantile(modelSkrim$depth_cm, c(0.25, 0.5, 0.75))

p_base <- ggplot(modelSkrim, aes(x = depth_cm)) +
  geom_histogram(binwidth = binwidth)
hist_data <- layer_data(p_base)
x_values <- c(mean_depth, quantiles_depth)
bin_indices <- findInterval(x_values, hist_data$xmin, rightmost.closed = TRUE)
y_values <- hist_data$count[bin_indices]

label_data <- data.frame(
  x = x_values,
  y = y_values,
  label = c(paste("mean:", round(mean_depth), "cm"),
            paste("Q1:", round(quantiles_depth[1]), "cm"),
            paste("median:", round(quantiles_depth[2]), "cm"),
            paste("Q3:", round(quantiles_depth[3]), "cm")
  ))

densitySkrim <- ggplot(modelSkrim, aes(x = depth_cm)) +
  geom_histogram(aes(fill = after_stat(x)), 
                 binwidth = binwidth, color = "black", size = 0.1) +
  scale_fill_viridis_c(limits = color_limits, direction = -1, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(modelSkrim$depth_cm), by = 200),
                     labels = scales::label_number(suffix = " cm")) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text_repel(data = label_data, aes(x = x, y = y, label = label),
                  direction = "y",
                  min.segment.length = 0,
                  nudge_x = max(hist_data$x) * 0.4,
                  nudge_y = max(hist_data$count) * 1,
                  size = 3) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.line.x = element_line()
  )

insetwidthOrskog <- 0.89 - 0.47 
scalingfactor <- max(modelSkrim$depth_cm) / max(modelOrskog$depth_cm)
pSkrim <- mapSkrim + labs(tag = "(a)") +
  inset_element(densitySkrim, 
                left = 0.12, # Not aligning to panel, seems a bug
                bottom = 0.03, 
                right = 0.12 + insetwidthOrskog * scalingfactor, 
                top = 0.5, 
                align_to = "full", clip = TRUE, ignore_tag = TRUE)

## Combined ####

pAll <- pSkrim + pOrskog + plot_layout(ncol = 1, nrow = 2, heights = c(2, 3))
ggsave(pAll, filename = 'map-distribution.pdf', path = "ms/figures",
       width =210-30, height = 297-30-40, units = 'mm') #A4 page 210x297

# Model metrics by faceting ####

library(tidyverse)

orskog <- read_csv("output/modelmetrics.csv") 
skrim <- read_csv("output/Skrim/modelmetrics.csv")
plotting <- bind_rows(orskog = orskog, skrim = skrim, .id = 'site') %>% 
  filter(.metric != 'rmse') %>% 
  mutate(
    site = fct_relevel(site, "orskog"),
    metric = case_when(
      .metric == "ccc" ~ "unitless\nConcordance correlation",
      .metric == "rsq" ~ "unitless\nR-squared",
      .metric == "mae" ~ "cm\nMean absolute error"),
    metric = fct_relevel(metric,
                         "unitless\nConcordance correlation",
                         "unitless\nR-squared",
                         "cm\nMean absolute error"),
    model = case_when(
      model == "DMK" ~ "DMK (2)",
      model == "Terrain" ~ "terrain (21)",
      model == "TerrainDMK" ~ "terrain + DMK (23)",
      model == "RadiometricTerrain" ~ "terrain + radiometric (25)",
      model == "RadiometricTerrainDMK" ~ "all predictors (27)"),
    model = fct_relevel(model,
                        "DMK (2)",
                        "terrain (21)",
                        "terrain + DMK (23)",
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
  facet_wrap(~metric, nrow = 1, scales = "free_x", axes ="margins",
             strip.position = "bottom") +
  scale_color_discrete(labels = c(orskog = "\u00D8rskogfjellet", 
                                  skrim = "Skrimfjella")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(-0.15, -0.17),
        legend.key.spacing.y = unit(-2, "mm"),
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave(g1, filename = 'modelmetrics-faceting.pdf', path = "ms/figures",
       width =210-30, height = (240-40)/2.5, units = 'mm') #copernicus.cls page 210x240

# Model metrics ####

library(tidyverse)
library(patchwork)

orskog <- read_csv("output/modelmetrics.csv") 
skrim <- read_csv("output/Skrim/modelmetrics.csv")
plotting <- bind_rows(orskog = orskog, skrim = skrim, .id = 'site') %>% 
  filter(.metric != 'rmse') %>% 
  mutate(
    site = fct_relevel(site, "skrim"),
    model = case_when(
      model == "Radiometric" ~ "radiometric (4)",
      model == "Terrain" ~ "terrain (21)",
      model == "DMK" ~ "DMK (2)",
      model == "RadiometricTerrain" ~ "terrain + radiometric (25)",
      model == "RadiometricDMK" ~ "radiometric + DMK (6)",
      model == "TerrainDMK" ~ "terrain + DMK (23)",
      model == "RadiometricTerrainDMK" ~ "all predictors (27)"),
    model = fct_relevel(model,
                        "DMK (2)",
                        "radiometric (4)",
                        "radiometric + DMK (6)",
                        "terrain (21)",
                        "terrain + DMK (23)",
                        "terrain + radiometric (25)",
                        "all predictors (27)"))
str(plotting)

cols <- c("orskog" = "#d95f02", "skrim" = "#1b9e77") #https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3

g1 <- plotting |> 
  filter(.metric == "ccc") |> 
  ggplot() +
  geom_linerange(aes(y = model, 
                     x = mean, 
                     xmin = mean - std_err, 
                     xmax = mean + std_err, 
                     color = site),
                 position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_point(aes(y = model, x = mean, color = site),
             position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_text(aes(y = model, x = mean, label = signif(mean, 2), group = site),
            position = position_dodge2(width = 1, reverse = TRUE),
            color = "grey50", size = 2.5) +
  scale_colour_manual(values = cols) +
  annotate(geom='text', x=0.05, y=6.85, size=3, label='\u00D8rskogfjellet', color = "#d95f02") +
  annotate(geom='text', x=0.05, y=7.15, size=3, label='Skrimfjella', color = "#1b9e77") +
  guides(color = "none") +
  labs(subtitle = "Concordance correlation") +
  xlab("unitless") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())

g2 <- plotting |> 
  filter(.metric == "rsq") |> 
  ggplot() +
  geom_linerange(aes(y = model, 
                     x = mean, 
                     xmin = mean - std_err, 
                     xmax = mean + std_err, 
                     color = site),
                 position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_point(aes(y = model, x = mean, color = site),
             position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_text(aes(y = model, x = mean, label = signif(mean, 2), group = site),
            position = position_dodge2(width = 1, reverse = TRUE),
            color = "grey50", size = 2.5) +
  scale_colour_manual(values = cols) +
  guides(col = "none") + 
  labs(subtitle = "R-squared") +
  xlab("unitless") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())

g3 <- plotting |> 
  filter(.metric == "mae") |> 
  ggplot() +
  geom_linerange(aes(y = model, 
                     x = mean, 
                     xmin = mean - std_err, 
                     xmax = mean + std_err, 
                     color = site),
                 position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_point(aes(y = model, x = mean, color = site),
             position = position_dodge2(width = 0.35, reverse = TRUE)) +
  geom_text(aes(y = model, x = mean, label = signif(mean, 2), group = site),
            position = position_dodge2(width = 1, reverse = TRUE),
            color = "grey50", size = 2.5) +
  scale_colour_manual(values = cols) +
  guides(color = "none") +
  labs(subtitle = "Mean absolute error") +
  xlab("cm") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())

g123 <- (g1 | g2 | g3) + plot_layout(axes = "collect_y")
ggsave(filename = 'modelmetrics.pdf', path = "ms/figures",
       width =210-30, height = (240-40)/2, units = 'mm') #copernicus.cls page 210x240
ggsave(filename = 'modelmetrics.png', path = "ms/figures",
       width =210-30, height = (240-40)/2, units = 'mm', dpi = 300) #copernicus.cls page 210x240

# Calibration plots ####

library(tidyverse)
library(patchwork)

orskog <- read_csv("output/calibrationplot.csv")
skrim <- read_csv("output/Skrim/calibrationplot.csv")

probs <- c(0.25, 0.5, 0.75)

scatter1 <- ggplot(orskog, aes(x = depth_cm, y = .pred)) +
  geom_point(size = 0.8, alpha=0.5) +
  geom_smooth(method = 'loess', formula= y ~ x, se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  tune::coord_obs_pred(ratio = 1) +
  labs(x = "Observed depth (cm)", y = "Predicted depth (cm)") +
  theme_bw()
dens <- density(orskog$depth_cm)
df <- data.frame(x = dens$x, y = dens$y)
quantiles <- quantile(orskog$depth_cm, prob=probs)
df$quant <- factor(findInterval(df$x, quantiles))
density1 <- ggplot(df, aes(x,y)) + 
  geom_area(aes(y = y, fill=quant)) +
  geom_line() +
  scale_fill_brewer(guide="none", palette = "Greys") +
  theme_void()

scatter2 <- ggplot(skrim, aes(x = depth_cm, y = .pred)) +
  geom_point(size = 0.8, alpha=0.5) +
  geom_smooth(method = 'loess', formula= y ~ x, se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  tune::coord_obs_pred(ratio = 1) +
  labs(x = "Observed depth (cm)", y = "Predicted depth (cm)") +
  theme_bw()
dens <- density(skrim$depth_cm)
df <- data.frame(x = dens$x, y = dens$y)
quantiles <- quantile(skrim$depth_cm, prob=probs)
df$quant <- factor(findInterval(df$x, quantiles))
density2 <- ggplot(df, aes(x,y)) + 
  geom_area(aes(y = y, fill=quant)) +
  geom_line() +
  scale_fill_brewer(guide="none", palette = "Greys") +
  theme_void()

density2 + density1 + scatter2 + scatter1 +
  plot_layout(nrow = 2, byrow = TRUE, heights = unit(c(-1, 4), c("null", "null"))) +
  plot_annotation(tag_levels =  list(c('(a)', '(b)')))
ggsave(filename = 'calibration_plots.pdf', path = "ms/figures",
       width =(210-30), height = (240-40)/2, units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'calibration_plots.png', path = "ms/figures",
#        width =(210-30), height = (240-40)/2, units = 'mm', dpi=300) #copernicus.cls page 210x240

# Variable importance ####

library(tidyverse)
library(patchwork)

orskog <- read_csv("output/variable_importance.csv") %>% 
  filter(type != "perm.ranger")
skrim <- read_csv("output/Skrim/variable_importance.csv") %>% 
  filter(type != "perm.ranger")

# Read correlation matrices for insets
corr_orskog <- read_csv("output/variable_importance_rank_correlations.csv")
corr_skrim <- read_csv("output/Skrim/variable_importance_rank_correlations.csv")

# Create simple correlation matrix plots for insets
create_simple_corr_plot <- function(corr_data) {
  corr_data %>%
    mutate(
      method1 = case_when(
        method1 == "firm" ~ "FIRM",
        method1 == "perm.vip" ~ "permutation", 
        method1 == "shap" ~ "Shapley",
        TRUE ~ method1
      ),
      method2 = case_when(
        method2 == "firm" ~ "FIRM",
        method2 == "perm.vip" ~ "permutation", 
        method2 == "shap" ~ "Shapley",
        TRUE ~ method2
      )
    ) %>%
    ggplot(aes(method1, method2)) +
    geom_tile(fill = "white", color = "black") +
    geom_text(aes(label = round(correlation, 2)), size = 2.5, color = "black") +
    coord_flip() +
    theme_void() +
    theme(axis.text = element_text(size = 6),
          plot.margin = margin(2, 2, 2, 2))
}

corr_plot_orskog <- create_simple_corr_plot(corr_orskog)
corr_plot_skrim <- create_simple_corr_plot(corr_skrim)

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
pOrskog <- plot.orskog %>% 
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
        axis.text.y.right = element_text(size = 8, color = "grey50", hjust = 0),
        legend.position = "inside",
        legend.position.inside = c(0.77, 0.22)) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(tag = "(b)") +
  inset_element(corr_plot_orskog, left = 0.7, bottom = 0.7, right = 1, top = 1, 
                align_to = "plot", clip = TRUE, ignore_tag = TRUE)
  
plot.skrim <- left_join(skrim, skrim.rank, by = join_by(Variable))
labs.skrim <- plot.skrim |> 
  distinct(Variable, removed, Imp.median) |>
  replace_na(list(removed = "")) |>
  arrange(Imp.median) |> 
  pull(removed) |> 
  stringr::str_wrap(width = 40)
pSkrim <- plot.skrim %>% 
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
        axis.text.y.right = element_text(size = 8, color = "grey50", hjust = 0)) +
  labs(tag = "(a)") +
  inset_element(corr_plot_skrim, left = 0.7, bottom = 0.1, right = 1, top = 0.4, 
                align_to = "plot", clip = TRUE, ignore_tag = TRUE)

## Combined ####

pAll <- pSkrim + pOrskog + plot_layout(ncol = 1, axes = 'collect') 

ggsave(pAll, filename = 'variable_importance.pdf', path = "ms/figures",
       width =(210-30), height = (240-40)/1.25, units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'variable_importance.png', path = "ms/figures",
#        width =(210-30), height = (240-40)/1.25, units = 'mm', dpi=300) #copernicus.cls page 210x240

# Partial dependence plots ####

library(tidyverse)
library(patchwork)

units <- read_csv("ms/tables/predictors.csv") |> 
  select(Code, Units) |> 
  bind_rows(tribble(~ Code, ~ Units,
                    "DMK_shallow", "binary",
                    "DMK_unknown", "binary"))

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
      TRUE ~ feature)) |> 
  left_join(units, by = c("feature" = "Code")) |> 
  mutate(
    feature = paste0(feature, " (", Units, ")"),
    feature = fct_reorder(feature, -Imp.median))

g1 <- ggplot(data = orskog, aes(x = .borders, y = .value, group = .id)) + 
  geom_line(data = filter(orskog, .type == 'ice'), alpha = 0.01) +
  geom_line(data = filter(orskog, .type == 'pdp'), color = "red") +
  geom_rug(data = filter(orskog, .type == "observed"), sides = "b") +
  coord_cartesian(ylim = c(0, 400)) +
  facet_wrap(~feature, scales = "free_x") + 
  ylab("Peat depth (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
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
      TRUE ~ feature)) |> 
  left_join(units, by = c("feature" = "Code")) |> 
  mutate(
    feature = paste0(feature, " (", Units, ")"),
    feature = fct_reorder(feature, -Imp.median))

g2 <- ggplot(data = skrim, aes(x = .borders, y = .value, group = .id)) + 
  geom_line(data = filter(skrim, .type == 'ice'), alpha = 0.01) +
  geom_line(data = filter(skrim, .type == 'pdp'), color = "red") +
  geom_rug(data = filter(skrim, .type == "observed"), sides = "b") +
  coord_cartesian(ylim = c(0, 250)) +
  facet_wrap(~feature, scales = "free_x") + 
  ylab("Peat depth (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())

g2 + g1 + plot_layout(ncol = 1, axes = 'keep') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'partial_dependence.pdf', path = "ms/figures",
       width =(210-30), height = (240-60), units = 'mm') #copernicus.cls page 210x240
# ggsave(filename = 'partial_dependence.png', path = "ms/figures",
#        width =(210-30), height = (240-60), units = 'mm', dpi = 300) #copernicus.cls page 210x240

# GPR wave velocity calibrations ####

library(tidyverse)
library(patchwork)

caldata <- read_csv("output/GPRcalibration-caldata.csv")
modelfit <- read_csv("output/GPRcalibration-modelfit.csv")

g1 <- caldata %>% 
  ggplot(aes(x=OWTT, y=depth_m)) +
  geom_point() +
  geom_ribbon(data=modelfit, aes(ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=modelfit) +
  geom_abline(slope = 0.03, intercept = 0, color="red", lty=2) +
  annotate("text", x = 20, y = 5, label = "R-squared = 0.947") +
  annotate("text", x = 100, y = 2.75, label = "0.03 m/ns (fresh water)", color="red", angle=21) +
  annotate("text", x = 110, y = 4.5, label = "0.043 m/ns", angle=29) +
  labs(x = "one-way travel time (ns)", y = "depth (m)") +
  coord_cartesian(expand = FALSE) +
  theme_bw()

caldata <- read_csv("output/Skrim/GPRcalibration-caldata.csv")
modelfit <- read_csv("output/Skrim/GPRcalibration-modelfit.csv")

g2 <- caldata %>% 
  ggplot(aes(x=OWTT, y=depth)) +
  geom_point(pch=16, size=2) +
  geom_ribbon(data=modelfit, aes(ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=modelfit) +
  geom_abline(slope = 0.03, intercept = 0, color="red", lty=2) +
  annotate("text", x = 15, y = 3.5, label = "R-squared = 0.874") +
  annotate("text", x = 60, y = 1.65, label = "0.03 m/ns (fresh water)", color="red", angle=20) +
  annotate("text", x = 60, y = 2.5, label = "0.039 m/ns", angle=26) +
  labs(x = "one-way travel time (ns)", y = "depth (m)") +
  coord_cartesian(expand = FALSE) +
  theme_bw()

(g2 + g1) + plot_layout(ncol = 1, axes = 'collect_x') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'GPRwavevelocity.pdf', path = "ms/figures",
       width =(210-30), height = ((240-30)), units = 'mm') #copernicus.cls page 210x240

# Semivariograms ####

library(tidyverse)
library(patchwork)
library(ggtext)

orskog <- read_csv("output/variogram-all-cutoff200.csv")
skrim <- read_csv("output/Skrim/variogram-all-cutoff200.csv")

g1 <- ggplot(orskog, aes(x = dist, y = gamma)) +
  geom_point() +
  labs(
    x = "Lag distance (m)",
    y = "Semivariance (cm^2^)"
  ) +
  theme_bw() +
  theme(axis.title.y = element_markdown())
g2 <- ggplot(skrim, aes(x = dist, y = gamma)) +
  geom_point() +
  labs(
    x = "Lag distance (m)",
    y = "Semivariance (cm^2^)"
  ) +
  theme_bw() +
  theme(axis.title.y = element_markdown())

g2 + g1 + plot_layout(ncol = 1, guides = 'collect', axes = 'collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') 
ggsave(filename = 'semivariograms.pdf', path = "ms/figures",
       width =(210-30), height = ((240-30)/2), units = 'mm') #copernicus.cls page 210x240

# sessionInfo ####

sessioninfo::session_info()
