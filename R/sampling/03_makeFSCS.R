# Script adapted from https://github.com/DickBrus/TutorialSampling4DSM
# Citation: Brus, D.J., 2019. Sampling for digital soil mapping: A tutorial supported by R scripts. Geoderma, 338, pp.464-480.

library(tidyverse)
library(sf)
library(terra)

#Read data with coordinates and other attributes
df <- read_csv(file="output/preparedSamplingMatrix.csv")
vars <- c('RAD_K','RAD_TC','RAD_Th','RAD_U','slope')
summary(df)
cor(select(df, all_of(vars)))
dfs <- df %>% 
  mutate(across(all_of(vars), ~ scale(.)[, 1]))
summary(dfs)

#Set number of sampling locations to be selected
N <- 160

#Compute clusters
#source("kmeanspp.R") does not solve convergence problem 'Quick-TRANSfer stage steps exceeded maximum'
set.seed(42)
myClusters <- kmeans(select(dfs, all_of(vars)), 
                     centers=N, iter.max=100, nstart=100)
if(myClusters$ifault == 4) {warning("kmeans failed to converge")}
write_rds(myClusters, "output/makeFSCS.rds")
# myClusters <- read_rds("output/makeFSCS.rds")

str(myClusters)
df$cluster <- myClusters$cluster
dfs$cluster <- myClusters$cluster

#Calculate reference distances
#Minimum distances to cluster centers
rdist.out <- fields::rdist(x1=myClusters$centers, x2=select(dfs, all_of(vars)))
ids.mindist <- apply(rdist.out, MARGIN=1, which.min)
mindist <- diag(rdist.out[,ids.mindist])

#Mean distances to cluster centers
meandist <- c()
for (i in seq_len(N)) {
  center <- myClusters$centers[i, ,drop=FALSE]
  cluster <- filter(dfs, cluster == i) %>% 
    select(all_of(vars)) %>% 
    as.matrix()
  meandist[i] <- mean(as.numeric(fields::rdist(x1=center, x2=cluster)))
}

#Restricting access by r.walk.cost
rwalkcost <- rast("data/rwalkcost.tif")
plot(rwalkcost < 1950)
maxcost <- 1950 #arbitrary scale

dfsa <- filter(dfs, rwalkcost <= maxcost)
rdist.out <- fields::rdist(x1=myClusters$centers, x2=select(dfsa, all_of(vars)))
ids.mindist <- apply(rdist.out, MARGIN=1, which.min)
mindista <- diag(rdist.out[,ids.mindist])
mySample <- dfsa[ids.mindist,c("x","y")] %>%
  left_join(df, by = join_by(x, y)) %>% 
  mutate(for.cluster = seq_len(N))
  
#Restricting access by r.walk.cost increases distances by %: 
#(using per cluster minimums and means as references)
(sum(mindista) - sum(mindist))/(sum(meandist) - sum(mindist))

plot(meandist, ylim=c(0,max(meandist)))
points(mindist, col='blue')
points(mindista, col='orange')

# Write sample to file
cov <- rast("output/samplingCovariates.tif")
sample <- st_as_sf(mySample, coords = c('x','y'), crs = crs(cov))

plot(cov['slope'], ext = ext(sample))
plot(st_geometry(sample), add=TRUE)

st_write(sample, "output/sampling.gpkg", layer = "FSCS", append=FALSE)

# Manual adjustments for access (not captured by r.walk.cost due to errors in land cover data etc.)
any(duplicated(pull(mySample, for.cluster)))

fc <- 160
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[11],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 85
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[7],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 137
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[2],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 26
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[3],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 146
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[3],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 9
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[4],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 125
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[3],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 100
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[2],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 65
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[2],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 67
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[2],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 94
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[4],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 58
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[2],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 90
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[4],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 148
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[9],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 7
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[11],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

fc <- 68
filter(mySample, for.cluster==fc)
dfsa[order(rdist.out[fc,])[1:20],]
replacement <- mutate(dfsa[order(rdist.out[fc,])[5],], for.cluster=fc)
mySample <- bind_rows(replacement, mySample)

# Remove original samples for those with replacements, and replace distance values with actual values
adjSample <- distinct(mySample, for.cluster, .keep_all = TRUE) %>% 
  arrange(for.cluster) %>% 
  select(x, y) %>% 
  left_join(df, by = join_by(x, y)) %>% 
  mutate(for.cluster = seq_len(N))

adjdfsa <- select(adjSample, x, y) %>% 
  left_join(dfsa, by = join_by(x, y))
adjrdist <- fields::rdist(x1=myClusters$centers, x2=select(adjdfsa, all_of(vars)))
distadj <- diag(adjrdist)

#Restricting access by r.walk.cost and manually increases distances by %: 
#(using per cluster minimums and means as references)
(sum(distadj) - sum(mindist))/(sum(meandist) - sum(mindist))
sum(distadj - mindist)/sum(meandist - mindist)
sum(distadj)/sum(mindist) #1.776308
sum(distadj)/sum(meandist) #0.4621858

plot(meandist, ylim=c(0,max(meandist)))
points(mindist, col='blue')
points(distadj, col='orange')

#Plot clusters and sampling points
ggplot(slice_sample(df, n=10000)) +
  geom_point(mapping=aes(y=RAD_TC, x=slope, colour=factor(cluster)), size = 0.5) +
  geom_point(data=adjSample, 
             mapping=aes(y=RAD_TC, x=slope, fill=factor(for.cluster)), size=3, pch=21) +
  scale_y_continuous(name = "RAD_TC") +
  scale_x_continuous(name = "slope") +
  theme(legend.position="none") +
  labs(title = "Selected cluster representative for sampling (n=160)",
       subtitle = "over 2D of 5D-covariate space within peatland land cover area")
ggsave(filename = "output/FSCSaccessible.svg", width = 15, height = 15, units = "cm")

#Write sample to file
adjsample <- st_as_sf(adjSample, coords = c('x','y'), crs = crs(cov))
st_write(adjsample, "output/sampling.gpkg", layer = "FSCSaccessible", append=FALSE)
