# After @maloneMethodsImproveUtility2019 and @sauretteDivergenceMetricsDetermining2023

library(tidyverse)
library(clhs)

#Read data with coordinates and other attributes
df <- read_csv(file="output/preparedSamplingMatrix.csv")
vars <- c('RAD_K','RAD_TC','RAD_Th','RAD_U','slope')
summary(df)
cor(select(df, all_of(vars)))

df %>% 
  select(all_of(vars)) %>% 
  map_int(nclass.FD)
nb <- 200

op <- par(mfrow = c(2,3))
for (i in vars) {
  pull(df, i) %>% 
    hist(breaks=nb, main = i)
}
par(op)

#?entropy::KL.empirical
# y1 = seq(1,20,1)
# y2 = c(seq(1,17,1),19,16,4)
# entropy::KL.empirical(y1, y2, unit="log2")
# philentropy::KL(rbind(y1/sum(y1),y2/sum(y2)), unit="log2")

dfcovs <- df[,vars]
# v1 <- dfcovs$RAD_K
# v2 <- s.dfcovs$RAD_K
kl.easy <- function(v1, v2, nb) {
  breaks <- seq(min(v1), max(v1), length.out = nb+1)
  y1 <- hist(v1, breaks = breaks, plot=FALSE)$counts
  y2 <- hist(v2, breaks = breaks, plot=FALSE)$counts
  philentropy::KL(rbind(y1, y2), unit="log2", est.prob = "empirical")
}

nseq <- c(seq(50,200,25), seq(250,500,50)) # cLHS sample size
its <- 10  # number internal iterations with each sample size number
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
write_csv(res, file="output/sampleSizeKL.csv")
# res <- read_csv("output/sampleSizeKL.csv")

ggplot(res) +
  geom_point(mapping = aes(x=n, y = KL))

fit <- nls(KL ~ SSasymp(n, yf, y0, log_alpha), data = res)
fit
res <- mutate(res, fitted = fitted(fit))

ggplot(res) +
  geom_point(mapping = aes(x=n, y = KL)) +
  geom_line(mapping = aes(x=n, y=fitted)) +
  xlim(0, 500) + ylim(0, 5)

#Apply fit
xx <- seq(1, 500, 1)
jj <- predict(fit, list(n=xx))
normalized = 1 - (jj-min(jj))/(max(jj)-min(jj))

x <- xx
y <- normalized

sn <- tibble(x, y) %>%
  filter(y >= 0.95) %>% 
  slice_head(n = 1) %>% 
  pull(x)
sn

plot(x, y, xlab="sample number", ylab= "normalised KL", type="l", lwd=2)
x1<- c(-1, 500); y1<- c(0.95, 0.95)
lines(x1,y1, lwd=2, col="red")

x2<- c(sn, sn); y2<- c(0, 1)
lines(x2,y2, lwd=2, col="red")

