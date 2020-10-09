library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(rmcorr)
source("funcoes.R")

################################################################################
## Correlacao entre limitação temporal e espacial
################################################################################

################################################################################
## Gráficos ##
################################################################################
##Medias e min-max (linhas cruzadas)  
ab.sp.rsl$ssl.meanj <- jitter(ab.sp.rsl$ssl.mean)
ab.sp.rsl$tsl.meanj <- jitter(ab.sp.rsl$tsl.mean)
plot(tsl.meanj ~ ssl.meanj, data = ab.sp.rsl, xlim=c(0,1), ylim=c(0,1))
with(ab.sp.rsl, segments(x0 = ssl.min, x1 = ssl.max, y0=tsl.meanj, y1=tsl.meanj))
with(ab.sp.rsl, segments(x0 = ssl.meanj, x1 = ssl.meanj, y0=tsl.min, y1=tsl.max))

## Todos os pontos coloridos por spp
abundants %>%
    ggplot(aes(ssl, tsl)) +
    geom_jitter(aes(color=species)) +
    scale_color_viridis_d() +
    theme_bw()

## Separado por ano
abundants %>%
    ggplot(aes(ssl, tsl)) +
    geom_jitter(aes(color=species)) +
    facet_grid(~year)

##Grafico tot40 X tot12 
plot(abundants2$ssl_tot40, abundants2$tsl_tot12, xlab="Temporal seed limitation", ylab="Spatial seed limitation")

################################################################################
## Calculos da correlação
################################################################################
## Correlação entre as ssl e tsl médio, levando em conta medidas repetidas
## Pacote rmcorr
sl.tl.rmcor <- rmcorr(participant = factor(species), measure1 = ssl, measure2 = tsl, dataset=abundants, CI="bootstrap")
print(sl.tl.rmcor)
plot(sl.tl.rmcor)

## Correlação entre as médias (talvez baste)
cor(ab.sp.rsl[,c("ssl.mean", "tsl.mean")])

##Correlações a cada ano
y1<-abundants%>%
    filter (year=="1")
cor(y1$ssl, y1$tsl)

y2<-abundants%>%
    filter (year=="2")
cor(y2$ssl, y2$tsl)

y3<-abundants%>%
    filter (year=="3")
cor(y3$ssl, y3$tsl)

##Correlação total (tot40 e tot12)
cor(tudo82spp$ssl_tot40, tudo82spp$tsl_tot12)

