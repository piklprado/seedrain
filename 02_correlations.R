library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(rmcorr)
source("funcoes.R")

################################################################################
## Correlation between temporal and spatial limitation (tsl x ssl)
################################################################################

################################################################################
## Plots ##
################################################################################
## Mean and ranges of tsl and ssl  
ab.sp.rsl$ssl.meanj <- jitter(ab.sp.rsl$ssl.mean)
ab.sp.rsl$tsl.meanj <- jitter(ab.sp.rsl$tsl.mean)
plot(tsl.meanj ~ ssl.meanj, data = ab.sp.rsl, xlim=c(0,1), ylim=c(0,1))
with(ab.sp.rsl, segments(x0 = ssl.min, x1 = ssl.max, y0=tsl.meanj, y1=tsl.meanj))
with(ab.sp.rsl, segments(x0 = ssl.meanj, x1 = ssl.meanj, y0=tsl.min, y1=tsl.max))

## Points colored by species
abundants %>%
    ggplot(aes(ssl, tsl)) +
    geom_jitter(aes(color=species)) +
    scale_color_viridis_d() +
    theme_bw()

## By year
abundants %>%
    ggplot(aes(ssl, tsl)) +
    geom_jitter(aes(color=species)) +
    facet_grid(~year)

##  ssl x tsl
plot(abundants2$ssl_tot40, abundants2$tsl_tot12,
     xlab="Temporal seed limitation",
     ylab="Spatial seed limitation")

################################################################################
## Correlations
################################################################################
## Between mean ssl and tsl, with repeated measurements
## Package rmcorr
sl.tl.rmcor <- rmcorr(participant = factor(species), measure1 = ssl, measure2 = tsl, dataset=abundants, CIs="bootstrap")
print(sl.tl.rmcor)
plot(sl.tl.rmcor)

## Pearson correlation between mean ssl and tsl
cor(ab.sp.rsl[,c("ssl.mean", "tsl.mean")])
cor.test(ab.sp.rsl$ssl.mean, ab.sp.rsl$tsl.mean)

## Pearson correlation each year
## Year 1
y1<-abundants%>%
    filter (year=="1")
cor(y1$ssl, y1$tsl)
cor.test(y1$ssl, y1$tsl)

## Year2
y2<-abundants%>%
    filter (year=="2")
cor(y2$ssl, y2$tsl)
cor.test(y2$ssl, y2$tsl)

## Year 3
y3<-abundants%>%
    filter (year=="3")
cor(y3$ssl, y3$tsl)
cor.test(y3$ssl, y3$tsl)

## Correlation for all species
cor(tudo82spp$ssl_tot40, tudo82spp$tsl_tot12)

