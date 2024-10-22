library(tidyverse)
library(sf)
library(terra)
library(gstat)
library(clhs)

# Study area

predictors <- rast("output/predictors.tif")
sa <- st_read("data/Orskogfjellet-site.gpkg", "mask_studyarea")
extract(predictors, sa, ID = FALSE) |> 
  summary()

# Spatial structure in peat depth ####

depthpts <- read_csv("data/depth/all.csv") |> 
  st_as_sf(coords = c('X', 'Y'), crs = 25832)

# From probe data only
depth_probe <- depthpts |> 
  filter(source == "probe") |> 
  as("Spatial")  
emp_variog <- variogram(depth_cm ~ 1, cutoff = 1000, data = depth_probe)
print(emp_variog)
plot(emp_variog)
emp_variog <- variogram(depth_cm ~ 1, cutoff = 200, data = depth_probe)
print(emp_variog)
plot(emp_variog)

# From GPR data only
depth_gpr <- depthpts |> 
  filter(source %in% c('gpr', 'wisen2021')) |> 
  as("Spatial")
emp_variog <- variogram(depth_cm ~ 1, cutoff = 200, width = 5, 
                        data = depth_gpr[sample(seq_len(nrow(depth_gpr)), 
                                                nrow(depth_gpr)*0.1),])
plot(emp_variog)
print(emp_variog) # Range of 75 m
sqrt(140*2) # Nugget at 2.5 m of 16 cm

# At modeling resolution: average depth in 10 m cells
depthcells <- st_read("output/modeling.gpkg", "dataframe")
depthcells_sp <- as(depthcells, "Spatial")
emp_variog <- variogram(depth_cm ~ 1, cutoff = 1000, width = 10, 
                        data = depthcells_sp)
plot(emp_variog)
print(emp_variog) # Range of 230 m
sqrt(9350*2) # Sill at 230 m of 136 cm
sqrt(760*2) # Nugget at 10 m of 40 cm
# Isolated but similar peatland basins in the landscape can explain semivariance decline beyond range. 

# Sample attributes ####

depthcells <- st_read("output/modeling.gpkg", "dataframe", stringsAsFactors = TRUE)
count(depthcells, ar5cover)
count(depthcells, ar5cover, ar5soil)
count(depthcells, dmkdepth)
depthcells |> 
  st_drop_geometry() |> 
  group_by(ar5cover) |> 
  summarize(n = n(), depth_cm = mean(depth_cm))
depthcells |> 
  st_drop_geometry() |> 
  filter(ar5cover %in% c(30,50,60)) |> 
  group_by(ar5cover, ar5soil) |> 
  summarize(n = n(), depth_cm = mean(depth_cm))
depthcells |> 
  st_drop_geometry() |> 
  filter(ar5cover == 60) |> 
  group_by(dmkdepth) |> 
  summarize(n = n(), depth_cm = mean(depth_cm))
  
# Representativeness of the sample ####
# Representativeness cf. mapped peatland in the study area

predictors <- rast("output/predictors.tif")
plot(predictors)

ar5.myr <- st_read("data/Orskogfjellet-site.gpkg", layer="mask_ar5myr")
st_crs(ar5.myr) == st_crs(predictors)
plot(predictors[['elevation']])
plot(ar5.myr, add=TRUE)

st_intersects(depthcells, ar5.myr, sparse=FALSE) |> 
  apply(1, any) |> 
  (\(x) !x)() |> 
  sum() # 438 depth cells centers that do not intersect AR5 myr

ar5cells <- mask(predictors, ar5.myr)
plot(ar5cells)
vars <- names(predictors)
df <- as.data.frame(ar5cells, na.rm=TRUE)
summary(df)

df %>% 
  select(all_of(vars)) %>% 
  map_int(nclass.FD)
nb <- 200

for (i in vars) {
  pull(df, i) %>% 
    hist(breaks=nb, main = i)
}

dfcovs <- df[,vars]
kl.easy <- function(v1, v2, nb) {
  breaks <- seq(min(v1), max(v1), length.out = nb+1)
  y1 <- hist(v1, breaks = breaks, plot=FALSE)$counts
  y2 <- hist(v2, breaks = breaks, plot=FALSE)$counts
  philentropy::KL(rbind(y1, y2), unit="log2", est.prob = "empirical")
}

nseq <- 10*2^c(0:8) # cLHS sample size
its <- 3  # number internal iterations with each sample size number
res <- tibble(n=rep(nseq, each=its), 
              replicate=rep(seq_len(its), length(nseq)))

calcKL <- function(n, replicate) {
  ss <- clhs(dfcovs, size = n, progress = T)
  s.dfcovs <- dfcovs[ss,]
  KL <- mean(map2_dbl(dfcovs, s.dfcovs, kl.easy, nb = nb))
}

set.seed(42)
KL <- pmap_dbl(res, calcKL, .progress = TRUE)
res <- res %>% 
  mutate(KL = KL)
write_csv(res, file="output/sampleSizeKL-25preds.csv")
# res <- read_csv("output/sampleSizeKL-25preds.csv")

ggplot(res) +
  geom_point(mapping = aes(x=n, y = KL))

fit <- scam::scam(KL ~ s(n, k=20, bs='mpd'), data = res)
fit
res <- mutate(res, fitted = fitted(fit))

ggplot(res) +
  geom_point(mapping = aes(x=n, y = KL)) +
  geom_line(mapping = aes(x=n, y=fitted)) 

#Apply fit
xx <- seq(range(nseq)[1],range(nseq)[2], length.out = 1000)
jj <- predict(fit, list(n=xx))
normalized = 1 - (jj-min(jj))/(max(jj)-min(jj))

x <- xx
y <- normalized

sn <- tibble(x, y) %>%
  filter(y >= 0.95) %>% 
  slice_head(n = 1) %>% 
  pull(x)
sn # 285

plot(x, y, xlab="cLHS sample size", ylab= "normalised KL", type="l", lwd=2)
x1<- c(-1, 4000); y1<- c(0.95, 0.95)
lines(x1,y1, lwd=2, col="red")

x2<- c(sn, sn); y2<- c(0, 1)
lines(x2,y2, lwd=2, col="red")

dfrealized <- extract(ar5cells, depthcells, ID = FALSE) |> 
  drop_na()
KLrealized <- mean(map2_dbl(dfcovs, dfrealized, kl.easy, nb = nb))
KLrealizednormalized = 1 - (KLrealized-min(jj))/(max(jj)-min(jj))

sn_realized <- tibble(x, y) %>%
  filter(y >= KLrealizednormalized) %>% 
  slice_head(n = 1) %>% 
  pull(x)
points(sn_realized, KLrealizednormalized) # c(212, 0.941)
nrow(dfrealized) # 1550

# Collinearity ####

df <- st_read("output/modeling.gpkg", "dataframe") |> 
  st_drop_geometry() |> 
  as_tibble()
metadata <- c("sourceProbe", "sourceGPR", "sourceWisen2021", "sourceMyrarkivet")
auxiliary <- c("ar5cover", "ar5soil", "dmkdepth")
sum(!complete.cases(df))

svg("output/exploratoryDataAnalysis-corr.svg", width = 12, height = 12)
corrplot::corrplot(cor(select(df, !any_of(c(metadata, auxiliary)))), 
                   method = 'number', 
                   type = 'upper',
                   diag = TRUE,
                   order = 'hclust', 
                   number.cex = 1,
                   addCoefasPercent=TRUE)
dev.off()

select(df, !any_of(c(metadata, auxiliary))) %>% 
  cor() %>% 
  round(2)

dflong <- df |> 
  select(!any_of(c(metadata, auxiliary))) |> 
  pivot_longer(cols = -1, names_to = "covariate")

ggplot(dflong) + 
  geom_point(aes(x = value, y = depth_cm), pch=20) + 
  facet_wrap(facets = vars(covariate), scales = "free_x") +
  theme_minimal() 
ggsave("output/exploratoryDataAnalysis-pairs.svg", 
       height = 10, width = 20, unit="in")
