library(tidyverse)
library(terra)
library(sf)
library(mlr3verse)

training <- st_read("output/modeling.gpkg", layer="training")
celldepth <- training |> 
  select(depth_cm)

dmk <- list.files(path = "data/DMK", pattern = '*.shp', recursive = TRUE, 
                  full.names = TRUE) 
dmk <- map(dmk, st_read)
dmkmyr <- bind_rows(dmk) |> 
  select(myr) |>
  filter(myr != 0) |> 
  mutate(depth_class = case_when(
    myr >= 30 ~ "djup myr",
    myr < 30 ~ "grunn myr"
  ))
st_write(dmkmyr, "data/Orskogfjellet-site.gpkg", "dmkmyr", append=FALSE)

celldepth <- st_transform(celldepth, st_crs(dmkmyr))
df <- st_join(celldepth, dmkmyr)

# Qualitative evaluation ####
summary <- df |> 
  st_drop_geometry() |> 
  group_by(depth_class) |> 
  summarise(n = n(), meandepth = mean(depth_cm, na.rm=TRUE))
summary # NA from outside DMK peatland polygons
plotdat <- left_join(df, summary, by = join_by(depth_class)) |>
  filter(!is.na(depth_class)) |> 
  mutate(myaxis = paste0(depth_class, "\n", "n=", n),
         myaxis = fct(myaxis))

ggplot(plotdat, aes(x=myaxis, y=depth_cm)) + 
  geom_violin(width=1.2) +
  geom_boxplot(width=0.1, color = "grey") +
  geom_hline(yintercept = 100, color = "red", lty = 2) +
  xlab("DMK Torvdybde") + ylab("Dybde (cm)") +
  annotate("text", x = 1.5, y = 100, color = "red", size = 3,
           label = "DMK \ngrense")
ggsave("output/DMKqualitative.png", width = 12, height = 9, unit = "cm")

# Quantitative evaluation ####
# RMSE is highly dependent on point estimate for each of the two classes:
df <- df |> 
  mutate(depth_class_cm = case_when(
    myr >= 30 ~ 150,
    myr < 30 ~ 65,
    TRUE ~ NA
  ))
mse <- df |> 
  st_drop_geometry() |> 
  mutate(error = (depth_class_cm - depth_cm)^2) |> 
  summarize(mse = mean(error, na.rm=TRUE))
mse
sqrt(mse)

df <- df |> 
  mutate(depth_class_cm = case_when(
    myr >= 30 ~ 100,
    myr < 30 ~ 50,
    TRUE ~ NA
  ))
mse <- df |> 
  st_drop_geometry() |> 
  mutate(error = (depth_class_cm - depth_cm)^2) |> 
  summarize(mse = mean(error, na.rm=TRUE))
mse
sqrt(mse)

# Use random CV to evaluate calibrated DMK depths
dftsk <- df |>
  filter(!is.na(depth_class)) |> 
  st_drop_geometry() |> 
  select(depth_cm, depth_class)
tsk_depth <- as_task_regr(dftsk, target="depth_cm")

lrn_rf <- lrn("regr.ranger", num.trees = 1000)
lrn_rf$train(tsk_depth)
lrn_rf$predict(tsk_depth)$score(msr("regr.rmse"))

rcv1010 <- rsmp("repeated_cv", repeats = 10, folds = 10)
set.seed(123)
rr <- resample(tsk_depth, lrn_rf, rcv1010)
rr$aggregate(msr("regr.rmse"))

# Session info
sessionInfo()
