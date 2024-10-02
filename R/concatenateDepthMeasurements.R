library(tidyverse)
library(sf)

# Probe data ####

probe <- read_csv("data/depth/concatenatedRawProbe.csv")
probe <- st_as_sf(probe, coords = c('utm32E', 'utm32N'), crs = 25832)

## Cleaning probe data ####

filter(probe, is.na(depth_cm)) # 2 of 160 sampling cells have soil altered by construction
filter(probe, !is.na(comment)) |>
  arrange(comment) |> 
  print(n=160)

# Remove depths which have changed since LIDAR and radiometric surveys in 2015. 
# Do not remove depths that are affected by infrastructure already in 2015.
# Checked all points with comments: "infill", "fyllmasse", "gravemasser", "road"
probeclean <- probe |> 
  filter(!is.na(depth_cm)) |> 
  filter(!(from == "dgps2" & ptname == 276 & comment == "fyllmasse, ikke torv")) |> 
  filter(!(from == "dgps2" & ptname == 251 & comment == "road, fyllmasse paa veg"))

arrange(probeclean, desc(HRMS))

# GPR data ####

gprpicks <- st_read("data/depth/GPRpicks.gpkg", "picks-round2")
gprpicks <- select(gprpicks, file = field_1, trace =field_4, TWTT = field_6) |> 
  st_transform(st_crs(probe))

gprlines <- gprpicks %>% 
  group_by(file) |> 
  arrange(trace) |>
  select(file) |>
  summarize(file = first(file), do_union=FALSE) |> 
  st_cast("LINESTRING")
gprlines |> 
  mutate(length = st_length(gprlines)) |> 
  pull(length) |> 
  sum()

# Distance threshold to divide trace sequences into separate lines
threshold_distance <- 10
buffers_sf <- st_buffer(gprpicks, dist = threshold_distance/2) |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_sf() |> 
  rowid_to_column("group_id")
gprlines <- gprpicks |> 
  st_join(buffers_sf) |> 
  arrange(trace) |>
  group_by(file, group_id) |>
  summarize(file = first(file), group_id = first(group_id), 
            do_union=FALSE, .groups = "drop") |> 
  st_cast("LINESTRING")
gprlines |> 
  mutate(length = st_length(gprlines)) |> 
  pull(length) |> 
  sum()

## Convert GPR TWTT to depth ####

### Spatially join probe to GPR ####

st_crs(probeclean) == st_crs(gprpicks)

joinNearest <- function(x, y) { 
  nearest <- st_nearest_feature(x, y)
  ynearest <- y[nearest,]
  x <- mutate(x, xydistance = st_distance(x, ynearest, by_element = TRUE))
  bind_cols(x, st_drop_geometry(ynearest))
}

probeclean <- probeclean %>% 
  mutate(probeID = seq_len(nrow(probeclean)), .before = 1)
gprpicks <- group_by(gprpicks, file)

# Closest trace to probe for each profile (one-to-many)
sub2m <- gprpicks %>% 
  group_split() %>% 
  map(function(a) joinNearest(x = probeclean, y = a)) %>% 
  bind_rows() %>% 
  mutate(xydistance = units::drop_units(xydistance)) %>% 
  filter(xydistance < 2)

nrow(sub2m) == nrow(distinct(sub2m, probeID, file))

sub2m %>% 
  group_by(probeID) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n=10)

caldata <- sub2m %>% 
  select(probeID, depth_cm, file, trace, TWTT, xydistance) |> 
  mutate(depth_m = depth_cm/100, OWTT = TWTT/2)

### Calibrate wave velocity ####
m0 <- lm(depth_m ~ OWTT + 0, data = caldata)
summary(m0) # OWTT  0.042726

m1 <- lm(depth_m ~ OWTT + 0, weights = 1/xydistance, data = caldata)
summary(m1) # OWTT  0.04199

modelsummary::modelsummary(list(m0, m1), gof_map = "r.squared")

modelfit <- data.frame(OWTT = seq(0, max(caldata$OWTT)*1.05, length.out = 100))
modelfit <- modelfit %>% 
  bind_cols(predict(m0, modelfit, se.fit = TRUE, interval = "confidence")$fit) %>% 
  rename(depth_m = fit)

calplot <- caldata %>% 
  ggplot(aes(x=OWTT, y=depth_m)) +
  geom_point(aes(alpha = xydistance), pch=16, size=2) +
  geom_ribbon(data=modelfit, aes(ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=modelfit) +
  geom_abline(slope = 0.03, intercept = 0, color="red", lty=2) +
  scale_alpha(name = "probe-GPR lateral separation", range = c(1, 0.1), labels = scales::label_number(suffix = "m")) +
  annotate("text", x = 100, y = 1, label = "R-squared = 0.947") +
  annotate("text", x = 100, y = 2.75, label = "0.03 m/ns (fresh water)", color="red", angle=23) +
  annotate("text", x = 110, y = 4.25, label = "0.043 m/ns", angle=29) +
  labs(x = "one-way travel time (ns)", y = "depth (m)") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  theme(legend.position = c(0.3, 0.8),
        legend.background = element_blank())

ggsave('output/calibrationGPR.svg', calplot, width = 119, height = 84, units = "mm")
v <- units::set_units(m0$coefficients[[1]], "m/ns")

### TWTT to depth ####
gpr <- gprpicks |> 
  mutate(depth_cm = TWTT/2*units::drop_units(v)*100)

# Wisen 2021 (Impakt Geofysikk) ####

wisen1 <- read_tsv("data/depth/Wisen2021/Pick0_Torv_Myr_NTM06_2020_2021_combined.txt") |> 
  select(line = Line, trace =Trace, x_NTM06=x_NTM06, y_NTM06=y_NTM06, depth = `Depth @0.04`)
wisen2 <- read_csv("data/depth/Wisen2021/Pick0_Torv_Myr_02_NTM6.txt") |> 
  select(line = Name, trace = Trace, x_NTM06 = X_NTM6, y_NTM06 = Y_NTM6, depth = Depth)
wisen2021 <- list(wisen1, wisen2) |>
  set_names(nm = c("Pick0_Torv_Myr_NTM06_2020_2021_combined", "Pick0_Torv_Myr_02_NTM6")) |> 
  bind_rows(.id = "file") |> 
  mutate(depth_cm = depth*100) |> 
  select(-depth) |> 
  st_as_sf(coords = c('x_NTM06','y_NTM06'), crs=5106)

dupcheck <- st_equals(wisen2021)
wisen2021 <- mutate(wisen2021, hasSpatDuplicate = map_lgl(dupcheck, \(x) length(x) > 1))

# Filter out duplicates, because some show discrepancies suggestive of error (@wisenE39Romsdalsfjorden202020212021)
wisen2021 <- filter(wisen2021, !hasSpatDuplicate)

# Myrarkivet ####

myrarkivet <- st_read("data/depth/Myrarkivet.gpkg", "MR175-borekart1984")

# Concatenate sources ####

all <- list(probe = select(probeclean, depth_cm),
            gpr = select(gpr, depth_cm),
            wisen2021 = select(wisen2021, depth_cm),
            myrarkivet = select(myrarkivet, depth_cm))
all <- all |> 
  map(\(x) st_transform(x, crs = 25832)) |>
  map(\(x) st_sf(st_set_geometry(x, NULL), geometry = st_geometry(x))) |> 
  bind_rows(.id="source")

all |> 
  st_drop_geometry() |> 
  group_by(source) |> 
  summarize(n = n(), meandepth=mean(depth_cm))

st_delete("data/depth/all.csv")
st_write(all, "data/depth/all.csv", layer_options = "GEOMETRY=AS_XY", append=FALSE)
