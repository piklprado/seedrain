library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
source("funcoes.R")

################################################################################
## Data preparation ##
################################################################################

################################################################################
## 1. 49 species identified up to species level ##
################################################################################
raw <- read.csv("data/table_final_traps_months_traits.csv", as.is=TRUE)
## covert sampling year in factor
raw$year<-as.factor(raw$year)
## factor for zoocoric dispersion
raw$zooc <- raw$syndr=="Zoocoria"

################################################################################
## Object "abundants": 
## The 31 most abundant species ( species with at least 5 seeds sampled)
################################################################################
abundants <- raw[raw$most_abund,]
## Standard deviation of species variables (across species, n=31)
sd.mass<-sd(abundants$mass[abundants$year=="1"])
sd.log.mass<-sd(log(abundants$mass[abundants$year=="1"]))
sd.freq<-sd(abundants$freq_ad[abundants$year=="1"])
sd.Hloc<-sd(abundants$Hloc[abundants$year=="1"])
## Mean and medians too (convenience)
## Mean
mean.mass<-mean(abundants$mass[abundants$year=="1"])
mean.log.mass<-mean(log(abundants$mass[abundants$year=="1"]))
mean.freq<-mean(abundants$freq_ad[abundants$year=="1"])
mean.Hloc<-mean(abundants$Hloc[abundants$year=="1"])
## Median
median.mass<-median(abundants$mass[abundants$year=="1"])
median.log.mass<-median(log(abundants$mass[abundants$year=="1"]))
median.freq<-median(abundants$freq_ad[abundants$year=="1"])
median.Hloc<-median(abundants$Hloc[abundants$year=="1"])

## Variable standardization in the 'abundants" dataframe. 
abundants %<>%
    mutate(sp.i = gsub("\\W*\\b(\\w)\\w*?\\b\\W*", "\\1", species),
           log.mass = log(mass),
           mass2 = (mass - mean.mass)/sd.mass,
           log.mass2 = (log.mass - mean.log.mass)/sd.log.mass,
           freq2 = (freq_ad - mean.freq)/sd.freq,
           height2 = (Hloc - mean.Hloc)/sd.Hloc)

## Spatial limitation : proportion of seed collectors in which the species was not recorded.
abundants$ssl <- with(abundants, nt_a/nt_tot)

## Temporal limitation: sampling months at which the species was not recorded
abundants$tsl <- with(abundants, nm_a/nm_tot)

################################################################################
## Ranking of most abundant species, according mean spatial and temporal limitation
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

## Adds the ssl ranking to data , using seed mass to order ties
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$ssl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.ssl <- 1:nrow(ab.sp.rsl)

## Adds the tsl ranking to data , using seed mass to order ties
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$tsl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.tsl <- 1:nrow(ab.sp.rsl)

################################################################################
## 2. A data frame with all species and sampling years in separate columns
################################################################################

tudo82spp <- read.csv("data/all_82spp_years_in_columns_mean_total.csv", as.is=FALSE)
## Calculates tsl and ssl 
tudo82spp$ssl_tot40<-tudo82spp$nt_a_tot40/tudo82spp$nt_tot_tot40
tudo82spp$tsl_tot12<-tudo82spp$nm_a_tot12/tudo82spp$nm_tot_tot12
## A new dataframe with only the 31 most abundant species, from the object above
abundants2<-tudo82spp %>%
    filter(most_abund == TRUE) %>%
    mutate(sp.i = gsub("\\W*\\b(\\w)\\w*?\\b\\W*", "\\1", species))
## Variable standardization
abundants2$log.mass2 <- with(abundants2, (log(mass)-mean(log(mass)))/sd(log(mass)))
abundants2$freq2 <- with(abundants2, (freq_ad-mean(freq_ad))/sd(freq_ad))
abundants2$height2 <- with(abundants2, (Hloc-mean(Hloc))/sd(Hloc))

################################################################################
## Temporal distribution of seedrain
################################################################################
## Reading data
sinc <- read.csv("data/DadosBrutos_Assincronia.csv", as.is=TRUE)
## Some NA's that are actually zeroes
sinc[is.na(sinc)] <- 0
## Filtering species with at least one seed sampled in the collectors
sinc <- sinc[,c(TRUE, TRUE, apply(sinc[,-(1:2)],2,sum)>0)]
## A few other corrections in the data
sinc %<>%
    mutate(year_study = factor(year_study), month_study = factor(month_study))

## relative abundance of each species
## Only species with ate least 5 records
sp.abund <- gsub( " ", ".", unique(abundants$species))
sinc.ab <- sinc[, c(names(sinc)[1:2], sp.abund)]
## Abundances relative to the total abundances over the three years of sampling
sinc.ab.ar <- sweep(sinc.ab[, -(1:2)], MARGIN=2, STATS = apply(sinc.ab[, -(1:2)], 2, sum), FUN="/")
## Abundances relative to the total of seeds sampled each year
sinc.ab.ary <-
    sinc.ab %>%
    group_by(year_study) %>%
    mutate_at(3:ncol(sinc.ab), .funs = ab.rel) %>%
    as.data.frame()
