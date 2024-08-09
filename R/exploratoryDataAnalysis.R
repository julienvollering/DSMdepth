library(tidyverse)
library(sf)

# Spatial structure in peat depth ####

# Representativeness of the sample ####

# Collinearity ####

df <- st_read("output/modeling.gpkg", "dataframe") |> 
  st_drop_geometry() |> 
  as_tibble()
metadata <- c("sourceProbe", "sourceGPR", "sourceWisen2021", "sourceMyrarkivet")
any(!complete.cases(df))

select(df, !any_of(metadata)) |> 
  slice_sample(n = 1000) |> 
  pairs(pch='.')
corrplot::corrplot(cor(select(df, !any_of(metadata)) ), 
                   method = 'ellipse', diag = F)

dflong <- df |> 
  select(!any_of(metadata)) |> 
  pivot_longer(cols = -1, names_to = "covariate")

ggplot(dflong) + 
  geom_point(aes(x = value, y = depth_cm), pch=20) + 
  facet_wrap(facets = vars(covariate), scales = "free_x") +
  theme_minimal() 
ggsave("output/exploratoryDataAnalysis.svg", height = 10, width = 20, unit="in")
