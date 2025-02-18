library(tidyverse)
library(sf)

crssite <- st_crs("epsg:25833")

# Probe data ####

probe <- st_read("data/Skrim/Skrim-site.gpkg", "depth")
filter(probe, !is.na(note))
probe <- pivot_longer(probe, cols = starts_with("probe"), 
                      names_to = "replicate",
                      values_to = "depth_cm")
filter(probe, is.na(depth_cm)) 

probe |> 
  filter(grepl("plus|\\+", note))

st_write(probe, "data/Skrim/depth_probe.csv", 
         layer_options = "GEOMETRY=AS_XY", delete_layer = TRUE)

# GPR data ####

## Calibrate wave velocity ####
caldata <- read_csv("data/Skrim/GPRcalibration.csv")

m0 <- lm(depth ~ OWTT + 0, data = caldata)
summary(m0) # OWTT  0.03867

modelfit <- data.frame(OWTT = seq(0, max(caldata$OWTT)*1.05, length.out = 100))
modelfit <- modelfit %>% 
  bind_cols(predict(m0, modelfit, se.fit = TRUE, interval = "confidence")$fit) %>% 
  rename(depth = fit)
mean(abs(m0$residuals)) #MAE (m) 0.2941374

caldata %>% 
  ggplot(aes(x=OWTT, y=depth)) +
  geom_point(pch=16, size=2) +
  geom_ribbon(data=modelfit, aes(ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=modelfit) +
  geom_abline(slope = 0.03, intercept = 0, color="red", lty=2) +
  annotate("text", x = 60, y = 3.5, label = "R-squared = 0.874") +
  annotate("text", x = 60, y = 1.7, label = "0.03 m/ns (fresh water)", color="red", angle=22) +
  annotate("text", x = 60, y = 2.5, label = "0.039 m/ns", angle=28) +
  labs(x = "one-way travel time (ns)", y = "depth (m)") +
  coord_cartesian(expand = FALSE) +
  theme_bw()

write_csv(caldata, "output/Skrim/GPRcalibration-caldata.csv")
write_csv(modelfit, "output/Skrim/GPRcalibration-modelfit.csv")

v <- units::set_units(m0$coefficients[[1]], "m/ns")

## TWTT to depth ####
mire1 <- st_read("data/Skrim/GPRpicks.gpkg", "mire1")
mire2 <- st_read("data/Skrim/GPRpicks.gpkg", "mire2")
mire3 <- st_read("data/Skrim/GPRpicks.gpkg", "mire3")
gpr <- bind_rows(mire1, mire2, mire3, .id = "mire") |> 
  mutate(depth_cm = TWT/2*units::drop_units(v)*100) |> 
  select(mire, Name, depth_cm)

gprlines <- bind_rows(mire1, mire2, mire3) %>% 
  group_by(Name) |> 
  arrange(Trace) |>
  select(Name) |>
  summarize(Name = first(Name), do_union=FALSE) |> 
  st_cast("LINESTRING")
gprlines |> 
  mutate(length = st_length(gprlines)) |> 
  pull(length) |> 
  sum()

# Distance threshold to divide trace sequences into separate lines
threshold_distance <- 10
buffers_sf <- st_buffer(gpr, dist = threshold_distance/2) |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_sf() |> 
  rowid_to_column("group_id")
gprlines <- gpr |> 
  st_join(buffers_sf) |> 
  group_by(Name, group_id) |>
  summarize(do_union=FALSE, .groups = "drop") |> 
  st_cast("LINESTRING")
gprlines |> 
  mutate(length = st_length(gprlines)) |> 
  pull(length) |> 
  sum()

# Concatenate sources ####

all <- list(probe = select(probe, depth_cm),
            gpr = select(gpr, depth_cm))
all <- all |> 
  map(\(x) st_transform(x, crs = crssite)) |>
  map(\(x) st_sf(st_set_geometry(x, NULL), geometry = st_geometry(x))) |> 
  bind_rows(.id="source")

all |> 
  st_drop_geometry() |> 
  group_by(source) |> 
  summarize(n = n(), meandepth=mean(depth_cm))

st_delete("data/Skrim/depth_all.csv")
st_write(all, "data/Skrim/depth_all.csv", layer_options = "GEOMETRY=AS_XY", append=FALSE)
