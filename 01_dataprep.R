library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
source("funcoes.R")

################################################################################
## Preparacao dos dados ##
################################################################################

################################################################################
## 1. Planilha somente com as 49 spp identificadas até o nível de spp ##
################################################################################
raw <- read.csv("data/table_final_traps_months_traits.csv", as.is=TRUE)
##transformando ano em fator no objeto raw
raw$year<-as.factor(raw$year)
## Variavel de Zoocoria no objeto raw
raw$zooc <- raw$syndr=="Zoocoria"

################################################################################
## Objeto "abundants": 
## Selecionando apenas as 31 especies mais abundantes no objeto raw: pelo menos 5 registros
################################################################################
abundants <- raw[raw$most_abund,]

##Criando um objeto sd usando apenas um dos anos (os valores para os três anos são iguais, então
## escolhi o ano 1). Fizemos isso para que o "n" usado no cálculo de sd seja o número de espécies (aqui = 31)
## e não o número de linhas (que é 93 nesse arquivo, pois cada espécie tem o valor repetido para os três anos). 
sd.mass<-sd(abundants$mass[abundants$year=="1"])
sd.log.mass<-sd(log(abundants$mass[abundants$year=="1"]))
sd.freq<-sd(abundants$freq_ad[abundants$year=="1"])
sd.Hloc<-sd(abundants$Hloc[abundants$year=="1"])
## medias e medianas tambem. Não precisava, apenas conveniencia
mean.mass<-mean(abundants$mass[abundants$year=="1"])
mean.log.mass<-mean(log(abundants$mass[abundants$year=="1"]))
mean.freq<-mean(abundants$freq_ad[abundants$year=="1"])
mean.Hloc<-mean(abundants$Hloc[abundants$year=="1"])
## Medianas
median.mass<-median(abundants$mass[abundants$year=="1"])
median.log.mass<-median(log(abundants$mass[abundants$year=="1"]))
median.freq<-median(abundants$freq_ad[abundants$year=="1"])
median.Hloc<-median(abundants$Hloc[abundants$year=="1"])

## Variaveis padronizadas no objeto abundants
abundants %<>%
    mutate(sp.i = gsub("\\W*\\b(\\w)\\w*?\\b\\W*", "\\1", species),
           log.mass = log(mass),
           mass2 = (mass - mean.mass)/sd.mass,
           log.mass2 = (log.mass - mean.log.mass)/sd.log.mass,
           freq2 = (freq_ad - mean.freq)/sd.freq,
           height2 = (Hloc - mean.Hloc)/sd.Hloc)

## Calcula limitacao espacial (proporcao de coletores em que nao esteve) no objeto abundants
abundants$ssl <- with(abundants, nt_a/nt_tot)

## Calcula a limitacao temporal (meses do ano em que não esteve) no objeto abundants
abundants$tsl <- with(abundants, nm_a/nm_tot)

################################################################################
## Ranking das especies mais abundantes, de acordo com o ssl medio
################################################################################
ab.sp.rsl <- abundants %>%
    group_by(species) %>%
    summarise(ssl.mean = mean(ssl), ssl.sd = sd(ssl), ssl.min=min(ssl), ssl.max = max(ssl), ssl.tot = sum(ssl),
              tsl.mean = mean(tsl), tsl.sd = sd(tsl), tsl.min=min(tsl), tsl.max = max(tsl), tsl.tot = sum(tsl),
              l.ssl.mean = log(ssl.mean / (1-ssl.mean)), l.ssl.sd = log(ssl.sd / (1-ssl.sd)),
              l.ssl.lower = l.ssl.mean - l.ssl.sd, l.ssl.upper = l.ssl.mean + l.ssl.sd,
              l.tsl.mean = log(tsl.mean / (1-tsl.mean)), l.tsl.sd = log(tsl.sd / (1-tsl.sd)),
              l.tsl.lower = l.tsl.mean - l.tsl.sd, l.tsl.upper = l.tsl.mean + l.tsl.sd,            
              height = mean(Hloc), freq = mean(freq_ad), mass = mean(mass)) %>%
    ungroup() %>%
    mutate(sp.i = gsub("\\W*\\b(\\w)\\w*?\\b\\W*", "\\1", species))

## Adiciona ranking de ssl usando massa como desempate
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$ssl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.ssl <- 1:nrow(ab.sp.rsl)

## Adiciona ranking de tsl usando massa como desempate
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$tsl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.tsl <- 1:nrow(ab.sp.rsl)

################################################################################
## 2. Planilha com todas as espécies e com cada ano em colunas diferentes
################################################################################

tudo82spp <- read.csv("data/all_82spp_years_in_columns_mean_total.csv", as.is=FALSE)
##Criando as variáveis ssl e tsl para tot40 e tot12, respectivamente
tudo82spp$ssl_tot40<-tudo82spp$nt_a_tot40/tudo82spp$nt_tot_tot40
tudo82spp$tsl_tot12<-tudo82spp$nm_a_tot12/tudo82spp$nm_tot_tot12
##Criando um arquivo com somente as 31 espécies mais abundantes, a partir do arquivo tudo82spp
abundants2<-tudo82spp %>%
    filter(most_abund == TRUE) %>%
    mutate(sp.i = gsub("\\W*\\b(\\w)\\w*?\\b\\W*", "\\1", species))
## Padronizando as variáveis para esse novo objeto
abundants2$log.mass2 <- with(abundants2, (log(mass)-mean(log(mass)))/sd(log(mass)))
abundants2$freq2 <- with(abundants2, (freq_ad-mean(freq_ad))/sd(freq_ad))
abundants2$height2 <- with(abundants2, (Hloc-mean(Hloc))/sd(Hloc))

################################################################################
## Temporal distribution of seedrain
################################################################################
## Leitura dos dados
sinc <- read.csv("data/DadosBrutos_Assincronia.csv", as.is=TRUE)
## Ai caramba, onde devia ser zero tá como NA
sinc[is.na(sinc)] <- 0
## Ai caramba, tem especie com zero sementes
sinc <- sinc[,c(TRUE, TRUE, apply(sinc[,-(1:2)],2,sum)>0)]
## Outras correcoes
sinc %<>%
    mutate(year_study = factor(year_study), month_study = factor(month_study))

## Dados em abundância relativa de cada espécie 
## Apenas especies com pelo menos 5 sementes
sp.abund <- gsub( " ", ".", unique(abundants$species))
sinc.ab <- sinc[, c(names(sinc)[1:2], sp.abund)]
## Abundancias relativas supra-anual: em relação ao total dos 3 anos
sinc.ab.ar <- sweep(sinc.ab[, -(1:2)], MARGIN=2, STATS = apply(sinc.ab[, -(1:2)], 2, sum), FUN="/")
## Abundancias relativas ao total em cada ano
sinc.ab.ary <-
    sinc.ab %>%
    group_by(year_study) %>%
    mutate_at(3:ncol(sinc.ab), .funs = ab.rel) %>%
    as.data.frame()
