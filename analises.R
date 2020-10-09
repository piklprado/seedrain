library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(lme4)
library(bbmle)
library(rptR)
library(MuMIn)
library(rmcorr)
library(effects)
library(merTools)
library(parallel)
library(car)
source("funcoes.R")

## Preparacao dos dados ##
##Planilha somente com as 49 spp identificadas até o nível de spp
raw <- read.csv("table_final_traps_months_traits.csv", as.is=TRUE)
summary(raw)

##Com todas as espécies e com cada ano em colunas diferentes
tudo82spp <- read.csv("all_82spp_years_in_columns_mean_total.csv", as.is=FALSE)

##Criando as variáveis ssl e tsl para tot40 e tot12, respectivamente
tudo82spp$ssl_tot40<-tudo82spp$nt_a_tot40/tudo82spp$nt_tot_tot40
tudo82spp$tsl_tot12<-tudo82spp$nm_a_tot12/tudo82spp$nm_tot_tot12

summary(tudo82spp)

##Criando um arquivo com somente as 31 espécies mais abundantes, a partir do arquivo tudo82spp
abundants2<-tudo82spp %>%
 filter(most_abund == TRUE)

##transformando ano em fator no objeto raw
raw$year<-as.factor(raw$year)

## Variavel de Zoocoria no objeto raw
raw$zooc <- raw$syndr=="Zoocoria"

## Selecionando apenas as 31 especies mais abundantes no objeto raw: pelo menos 5 registros
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
    mutate(log.mass = log(mass),
           mass2 = (mass - mean.mass)/sd.mass,
           log.mass2 = (log.mass - mean.log.mass)/sd.log.mass,
           freq2 = (freq_ad - mean.freq)/sd.freq,
           height2 = (Hloc - mean.Hloc)/sd.Hloc)

## Calcula limitacao espacial (proporcao de coletores em que nao esteve) no objeto abundants
abundants$ssl <- with(abundants, nt_a/nt_tot)

## Calcula a limitacao temporal (meses do ano em que não esteve) no objeto abundants
abundants$tsl <- with(abundants, nm_a/nm_tot)

## Ranking das especies de acordo com o ssl medio
ab.sp.rsl <- abundants %>%
    group_by(species) %>%
    summarise(ssl.mean = mean(ssl), ssl.sd = sd(ssl), ssl.min=min(ssl), ssl.max = max(ssl),
              tsl.mean = mean(tsl), tsl.sd = sd(tsl), tsl.min=min(tsl), tsl.max = max(tsl),
              l.ssl.mean = log(ssl.mean / (1-ssl.mean)), l.ssl.sd = log(ssl.sd / (1-ssl.sd)),
              l.ssl.lower = l.ssl.mean - l.ssl.sd, l.ssl.upper = l.ssl.mean + l.ssl.sd,
              l.tsl.mean = log(tsl.mean / (1-tsl.mean)), l.tsl.sd = log(tsl.sd / (1-tsl.sd)),
              l.tsl.lower = l.tsl.mean - l.tsl.sd, l.tsl.upper = l.tsl.mean + l.tsl.sd,            
              height = mean(Hloc), freq = mean(freq_ad), mass = mean(mass)) %>%
    ungroup()

## Adiciona ranking de ssl usando massa como desempate
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$ssl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.ssl <- 1:nrow(ab.sp.rsl)

## Adiciona ranking de tsl usando massa como desempate
ab.sp.rsl <- ab.sp.rsl[order(-ab.sp.rsl$tsl.mean, -ab.sp.rsl$mass),]
ab.sp.rsl$rank.tsl <- 1:nrow(ab.sp.rsl)

##Como ficaram os quatro objetos manipulados aqui
summary(raw)
summary(abundants)
summary(ab.sp.rsl)
summary(abundants2)

#####################################################################################################################
## Gráficos de distribuição do número bruto de coletores e meses em que as 82 spp espécies estiveram presentes
#####################################################################################################################

##Para número de coletores (traps) que cada espécie esteve presente
par(mfrow=c(1,4), bty="n", cex.lab=1.5, mar=c(5,5,1,0))
boxplot (tudo82spp$nt_p1, ylab="Number of traps", xlab="Year 1", bty = "l", col=NULL, width=2)
boxplot (tudo82spp$nt_p2, ylab=NULL, xlab="Year 2", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nt_p3, ylab=NULL, xlab="Year 3", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nt_p_mean, ylab=NULL, xlab="Mean", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
par(mfrow=c(1,1), bty="o")

##Para número de meses que cada espécie esteve presente
par(mfrow=c(1,4), bty="n", cex.lab=1.5, mar=c(5,5,1,0))
boxplot (tudo82spp$nm_p1, ylab="Number of months", xlab="Year 1", bty = "l", col=NULL, width=2)
boxplot (tudo82spp$nm_p2, ylab=NULL, xlab="Year 2", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nm_p3, ylab=NULL, xlab="Year 3", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nm_p_mean, ylab=NULL, xlab="Mean", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
par(mfrow=c(1,1), bty="o")

##Para número de coletores (traps) juntando os 3 anos (tot40)
boxplot (tudo82spp$nt_p_tot40, ylab="Number of traps", xlab="Total", bty = "l", col=NULL, width=2)

##Para número de meses juntando os 3 anos (tot12)
boxplot (tudo82spp$nm_p_tot12, ylab="Number of months", xlab="Total", bty = "l", col=NULL, width=2)


################################################################################
## Correlacao entre limitação temporal e espacial
################################################################################
## Gráficos ##
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


## Calculo da correlação
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

##Grafico tot40 X tot12 
plot(abundants2$ssl_tot40, abundants2$tsl_tot12, xlab="Temporal seed limitation", ylab="Spatial seed limitation")

################################################################################
## Limitação espacial
################################################################################
## Graficos exploratorios ##
## SSLs ordenados e com maximo e minimo
ab.sp.rsl %>%
    arrange(rank.ssl) %>%
    ggplot(aes(rank.ssl, ssl.mean)) +
    geom_point() +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max)) +
    scale_x_continuous(breaks=1:nrow(ab.sp.rsl), labels= ab.sp.rsl$species[order(ab.sp.rsl$rank.ssl)]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1))

## Relacoes com as variaveis preditoras    
## Medias e min-max de ssl
ab.sp.rsl %>%
    mutate(mass=jitter(mass)) %>%
    ggplot(aes(mass, ssl.mean, color=species)) +
    geom_point(size=2) +
    ## scale_x_log10() +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max))

ab.sp.rsl %>%
    mutate(height=jitter(height)) %>%
    ggplot(aes(height, ssl.mean, color=species)) +
    geom_point(size=2) +
    ## scale_x_log10() +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max))

ab.sp.rsl %>%
    mutate(freq=jitter(freq)) %>%
    ggplot(aes(freq, ssl.mean, color=species)) +
    geom_point(size=2) +
    ## scale_x_log10() +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max))

### Ajuste do modelo, modelo médio e tudo mais ##
## Modelo inicial para dredge
ab.sl.full <- glmer(cbind(nt_a, nt_p) ~ log.mass2 + freq2 + height2 +
                    log.mass2:freq2 + log.mass2:height2 +
                    freq2:height2 + (1|species),
                    family=binomial,
                 data=abundants,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

##Usando a função dredge
options(na.action = "na.fail")
ab.sl.full.d <- dredge(ab.sl.full, beta="none")

## Modelos com deltaAIC <2
subset(ab.sl.full.d, delta < 2) 

## modelos selecionados
ab.sl.selected <- get.models(ab.sl.full.d, delta < 2)

## ICs dos efeitos de cada modelo selecionado
ab.sl.selected.IC <- lapply(ab.sl.selected, confint)

## Resulta numa lista de ICs
ab.sl.selected.IC

## VIF: nada de muito sério
vif(ab.sl.selected[[3]])

############### Modelo medio ################
sl.mavg <- model.avg(ab.sl.selected)
summary(sl.mavg)

## CIs dos efeitos medios
confint(sl.mavg, full=TRUE) # para full
confint(sl.mavg, full=FALSE) # para conditional


## Calculo do intervalo de previsao com efeito fixo e aletorio, por
## bootstrap Primeiro criamos uma dataframe com os valores de altura a
## prever e uma espécie qualquer (não faz diferença qual, pq os
## efeitos aleatorios sao sorteados para cada observaçao) Os valores
## de height2 neste df são os mesmos usados no passo anterior para
## calcular os efeitos fixos
new.data.ssl <- expand.grid(
    freq2 = quantile(abundants$freq2, c(0.25, 0.75)) ,
    height2 = quantile(abundants$height2, c(0.25, 0.75)),
    log.mass2 = quantile(abundants$log.mass2, seq(0,1, by=0.05)),
    species = unique(abundants$species)[1]
)

## Funcao para repetir a cada simulaçao bootstrap: simplesmente o predict para cada valor de altura
f1 <- function(.) predict(., newdata = new.data.ssl)

## Super funcao que faz o boostrap e devolve uma lista cujo primerio elemento sao os previstos e seus ICs
## para cada combinacao das preditoras em newdata
ssl.bootMer <- ic.bootMer(lista.modelos = ab.sl.selected, newdata = new.data.ssl, nsim = 1000, parallel=TRUE, ncpus = 6)
ssl.pred <- ssl.bootMer$predicted[, -4]
## Converte valores padronizados em valores na escala original e logitos em probabilidades
## no dataframe de previstos
ssl.pred %<>%
    mutate(height = height2*sd.Hloc + mean.Hloc,
           freq = freq2*sd.freq + mean.freq,
           log.mass = log.mass2*sd.log.mass + mean.log.mass,
           mass = exp(log.mass),
           pfit = ilogit(fit),
           plower = ilogit(lower),
           pupper = ilogit(upper),
           height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq)),
           )

## O plot final: observados, previsto pelo modelo geral e IC dos fixos e total
## As especies estao divididas em 4 grupos, delimitados pelas medianas de frequencia e altura
## Os previstos são so calculados para a mediana de freq e altura em cada grupo
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq))) %>%
    ggplot(aes(mass, ssl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max, color=species)) +
    geom_line(aes(y = pfit), data = ssl.pred)+
    geom_ribbon(aes(y = pfit, ymin = plower, ymax = pupper), data = ssl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("SSL")  
p1
    
## Escala logito
p2 <- 
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq))) %>%
    ggplot(aes(mass, l.ssl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin= l.ssl.upper, ymax = l.ssl.lower, color=species)) +
    geom_line(aes(y = fit), data = ssl.pred)+
    geom_ribbon(aes(y = fit, ymin = lower, ymax = upper), data = ssl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("Logito SSL")
p2

## Apenas as linhas previstas, para avaliar o modelo
## Escala logito
p3 <- ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("SSL previsto (logito)")
p3

## Escala prob
p4 <- ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("SSL previsto")
p4
    
##para ver todos os graficos 
p1 ## barras são ssl minimo e máximo
p2 ## barras são 1 desvio-padrão do logito do desvio-padrão do ssl
p3
p4

######### Replicability (R-squared) para os modelos selecionados ##############
## Usando a funcao do MuMIn (acho que é o que precisa, pois dá o R2 condicional (fixed + aleatorios)
## Metodo delta bate com os valores usando o pacote do Nakagawa (abaixo)
## Modelo nulo, que será usado de referência
ab.sl.null <- glmer(cbind(nt_a, nt_p) ~ (1|species),
                    family=binomial,
                    data=abundants)
for(i in 1:length(ab.sl.selected))
    print(r.squaredGLMM(ab.sl.selected[[i]], null = ab.sl.null))


##PI: parei aqui 
###Dri: Não mexi nesses gráficos abaixo. Fui direto para os modelos, mais abaixo

################################################################################
## Limitação temporal
################################################################################
## Graficos exploratorios ##
## TSLs ordenados e com maximo e minimo
ab.sp.rsl %>%
    arrange(rank.tsl) %>%
    ggplot(aes(rank.tsl, tsl.mean)) +
    geom_point() +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max)) +
    scale_x_continuous(breaks=1:nrow(ab.sp.rsl), labels=ab.sp.rsl$species[order(ab.sp.rsl$rank.tsl)]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1) )

## Relacoes com as variaveis preditoras    
## Medias e min-max de tsl
ab.sp.rsl %>%
    mutate(mass=jitter(mass)) %>%
    ggplot(aes(mass, tsl.mean, color=species)) +
    geom_point(size=2) +
    ##scale_x_log10() +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max))

ab.sp.rsl %>%
    mutate(height=jitter(height)) %>%
    ggplot(aes(height, tsl.mean, color=species)) +
    geom_point(size=2) +
    ## scale_x_log10() +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max))

ab.sp.rsl %>%
    mutate(freq=jitter(freq)) %>%
    ggplot(aes(freq, tsl.mean, color=species)) +
    geom_point(size=2) +
    ## scale_x_log10() +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max))


###
###Dri começou a modificar aqui
###

## Modelo inicial para dredge
ab.tl.full <- glmer(cbind(nm_a, nm_p) ~ log.mass2 + freq2 + height2 +
                    mass2:freq2 + mass2:height2 +
                    freq2:height2 + (1|species),
                    family=binomial,
                 data=abundants,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
                                      
## Usando a função dredge 
options(na.action = "na.fail")
ab.tl.full.d <- dredge(ab.tl.full, beta="none")

## Modelos com deltaAIC <2
subset(ab.tl.full.d, delta<2) 

## modelos selecionados
ab.tl.selected <- get.models(ab.tl.full.d, subset=delta<2)

## Modelo medio
tl.mavg <- model.avg(ab.tl.selected)

## CIs dos efeitos medios
confint(tl.mavg, full = TRUE) # para full model 
confint(tl.mavg, full = FALSE) # para conditional model


## Calculo do intervalo de previsao com efeito fixo e aletorio, por
## bootstrap Primeiro criamos uma dataframe com os valores de altura a
## prever e uma espécie qualquer (não faz diferença qual, pq os
## efeitos aleatorios sao sorteados para cada observaçao) Os valores
## de height2 neste df são os mesmos usados no passo anterior para
## calcular os efeitos fixos
new.data.tsl <- expand.grid(
    freq2 = quantile(abundants$freq2, c(0.25, 0.75)) ,
    height2 = quantile(abundants$height2, c(0.25, 0.75)),
    log.mass2 = quantile(abundants$log.mass2, seq(0,1, by=0.05)),
    species = unique(abundants$species)[1]
)

## Funcao para repetir a cada simulaçao bootstrap: simplesmente o predict para cada valor de altura
f1 <- function(.) predict(., newdata = new.data.tsl)

## Super funcao que faz o boostrap e devolve uma lista cujo primerio elemento sao os previstos e seus ICs
## para cada combinacao das preditoras em newdata
tsl.bootMer <- ic.bootMer(lista.modelos = ab.tl.selected, newdata = new.data.tsl, nsim = 1000, parallel=TRUE, ncpus = 6)
tsl.pred <- tsl.bootMer$predicted[, -4]

## Converte valores padronizados em valores na escala original e logitos em probabilidades
## no dataframe de previstos
tsl.pred %<>%
    mutate(height = height2*sd.Hloc + mean.Hloc,
           freq = freq2*sd.freq + mean.freq,
           log.mass = log.mass2*sd.log.mass + mean.log.mass,
           mass = exp(log.mass),
           pfit = ilogit(fit),
           plower = ilogit(lower),
           pupper = ilogit(upper),
           height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq)),
           )

## O plot final: observados, previsto pelo modelo geral e IC dos fixos e total
## As especies estao divididas em 4 grupos, delimitados pelas medianas de frequencia e altura
## Os previstos são so calculados para a mediana de freq e altura em cada grupo
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq))) %>%
    ggplot(aes(mass, tsl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max, color=species)) +
    geom_line(aes(y = pfit), data = tsl.pred)+
    geom_ribbon(aes(y = pfit, ymin = plower, ymax = pupper), data = tsl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("TSL")  
p1

## Escala logito
p2 <- 
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste("Height < ", median.Hloc), paste("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste("Freq < ", median.freq), paste("Freq > ", median.freq))) %>%
    ggplot(aes(mass, l.tsl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin= l.tsl.upper, ymax = l.tsl.lower, color=species)) +
    geom_line(aes(y = fit), data = tsl.pred)+
    geom_ribbon(aes(y = fit, ymin = lower, ymax = upper), data = tsl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("Logito TSL")
p2


## Apenas as linhas previstas, para avaliar o modelo
## Escala logito
p3 <- tsl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("TSL previsto (logito)")
p3

## Escala prob
p4 <- tsl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("TSL previsto")
p4

##para ver todos os graficos 
p1 ## barras são tsl minimo e máximo
p2 ## barras são 1 desvio-padrão do logito do desvio-padrão do ssl
p3
p4

######### Replicability (R-squared) para os modelos selecionados ##############
## Usando a funcao do MuMIn (acho que é o que precisa, pois dá o R2 condicional (fixed + aleatorios)
## Metodo delta bate com os valores usando o pacote do Nakagawa (abaixo)
## Modelo nulo, que será usado de referência
ab.tl.null <- glmer(cbind(nm_a, nm_p) ~ (1|species),
                    family=binomial,
                    data=abundants)
for(i in 1:length(ab.tl.selected))
    print(r.squaredGLMM(ab.tl.selected[[i]], null = ab.tl.null))



###
###Dri: Parei aqui
###


######################################################################################################
#Novo modelo considerando os total dos três anos, considerando 40 coletores (tot40) e 12 meses (tot12) 
######################################################################################################

##Vamos usar o objeto abundants2, que foi criado antes (a partir do tudo82spp) e seleciona somente as 31 spp mais abundantes
summary(abundants2)


## Padronizando as variáveis para esse novo objeto
abundants2$mass2 <- with(abundants2, (mass-mean(mass))/sd(mass))
abundants2$freq2 <- with(abundants2, (freq_ad-mean(freq_ad))/sd(freq_ad))
abundants2$height2 <- with(abundants2, (Hloc-mean(Hloc))/sd(Hloc))


########### Limitação Espacial ###############
#O modelo agora será um glm, pois as espécies são independentes, uma vez que não há repetição nos anos
#Modelo inicial para Dredge (COM ZOOCORIA AINDA)
ab.ssl.tot40.full <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ mass2 + freq2 + height2 + zooc +
                    mass2:freq2 + mass2:height2 + mass2:zooc +
                    freq2:height2 + freq2:zooc +
                    height2:zooc,
                    family=binomial,
                    data=abundants2)
                    
## PI: Acho complicado que um modelo com 10 termos ter sido selecionado, para um conjunto de dados com 31 observações.
## Fiquei avaliando alternativas que diminuam o n de coeficientes do modelo selecionado
## Meu primeiro candidato é variável de  zoocoria. São só 4 espécies zoocóricas, e este termo participa de
## 4 coeficientes
## E acho que massa pode ficar em escala log, pela grande assimetria de distribuição
## Detalhe: o glm pode usar as variaveis originais. Não muda nada, mas facilita interpretação
##ENTÃO:

#Modelo inicial para Dredge (SEM ZOOCORIA)
ab.ssl.tot40.full <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ log(mass) + freq_ad + Hloc +
                    log(mass):freq_ad + log(mass):Hloc + freq_ad:Hloc,
                    family="binomial",
                    data=abundants2)

## Quasibinomial
## Funcoes  auxiliares
## receita em https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf
dfun <- function(object) {
    with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}
x.quasibinom <- function(...) {
    res <- quasibinomial(...)
    res$aic <- binomial(...)$aic
    res
}
## Modelo ajustado de tal forma que o dredge funfe com quasibinomial
ab.ssl.tot40.full.q <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ log(mass) + freq_ad + Hloc +
                               log(mass):freq_ad + log(mass):Hloc + freq_ad:Hloc,
                           family="x.quasibinom",
                           data=abundants2)

## Dredge
options(na.action = "na.fail")
ab.ssl.tot40.full.d <- dredge(ab.ssl.tot40.full, beta="none", rank=AICc)
ab.ssl.tot40.full.q.d <- dredge(ab.ssl.tot40.full.q, beta="none", rank = "QAICc", chat = dfun(ab.ssl.tot40.ful))

## Modelos com deltaAIC <2 
subset(ab.ssl.tot40.full.d, delta < 2)

## quasibinomial
subset(ab.ssl.tot40.full.q.d, delta < 2) 

## modelos selecionados
ab.ssl.tot40.selected <- get.models(ab.ssl.tot40.full.d, delta < 2)

## ICs dos efeitos de cada modelo selecionado
ab.ssl.tot40.selected.IC <- lapply(ab.ssl.tot40.selected, confint)

## Resulta numa lista de ICs
ab.ssl.tot40.selected.IC

## VIF: 
ab.ssl.tot40.single <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ mass2 + freq2 + height2 + zooc,
                    family=binomial, data=abundants2)

vif(ab.ssl.tot40.single)



##PI: buscando graficos para mostrar o modelo
## Usando as medianas de log da massa e da altura para dividir as especies em 4 grupos, como
## Mediana das variaveis de massa e H
Hloc.median <- median(abundants2$Hloc)
log.mass.median  <-  median(log(abundants2$mass))
freq.median  <-  median(abundants2$freq_ad)
## a analise de efeitos acima fez
m1 <- ab.ssl.tot40.selected[[1]] ## melhor modelo
summary(m1) ## sobredisperso
dfun(m1)

## Calculos dos previstos
newdata <- expand.grid(freq_ad = seq(min(abundants2$freq_ad), max(abundants2$freq_ad), length =30),
                       Hloc = quantile(abundants2$Hloc, c(0.25, 0.75)),
                       mass = quantile(abundants2$mass, c(0.25, 0.75)))
previstos <- as.data.frame(predict(m1, newdata=newdata, se.fit=TRUE, type = "response"))
## Juntando os previstos com a variaveis resposta
newdata <- cbind(newdata,previstos)
newdata  <- mutate(newdata,
                   Hloc.class = ifelse(Hloc < Hloc.median, paste("H < ", Hloc.median), paste("H > ", Hloc.median)),
                   mass.class = ifelse(log(mass) < log.mass.median,
                                       paste("Mass < ", exp(log.mass.median)), paste("Mass > ", exp(log.mass.median))),
                   lower = fit - se.fit,
                   upper = fit + se.fit)
    
## Primeira tentativa de plot

## O grafico: ssl em funcao de frequencia, em quatro grupos de especies, separados pelas medianas de
## Hloc e log(massa). As linhas são os previstos pelo modelo para cada valor de frequencia e para os valores
## de Hloc e Log(mass) medianos de cada grupo.

abundants3 <- abundants2 %>%
    mutate(Hloc.class = ifelse(Hloc < Hloc.median, paste("H < ", Hloc.median), paste("H > ", Hloc.median)),
           mass.class = ifelse(log(mass) < log.mass.median, paste("Mass < ", exp(log.mass.median)),
                               paste("Mass > ", exp(log.mass.median))),
           freq.class =ifelse(freq_ad < freq.median, paste("Freq < ", freq.median), paste("Freq > ", freq.median)) ,
           ssl_tot40 = nt_a_tot40/40, 
           logito.ssl = log(nt_a_tot40/nt_p_tot40))

## Grafico com freq em X
    ggplot(abundants3, aes(freq_ad)) +
        geom_point(aes(y=ssl_tot40, color=zooc)) +
        ##geom_line(data = newdata, aes(freq_ad, fit)) +
        ##geom_ribbon(data = newdata, aes(freq_ad, ymin=lower, ymax =upper), alpha =0.2) +
        facet_grid(Hloc.class ~ mass.class) +
        theme_bw()
## Grafico com log massa em X
ggplot(abundants3, aes(mass)) +
        geom_point(aes(y=ssl_tot40, color=zooc)) +
        ##geom_line(data = newdata, aes(freq_ad, fit)) +
        ##geom_ribbon(data = newdata, aes(freq_ad, ymin=lower, ymax =upper), alpha =0.2) +
        facet_grid(Hloc.class ~ freq.class) +
        scale_x_log10() +
    theme_bw()
## idem com logito
## Grafico com log massa em X
ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=logito.ssl, color=zooc)) +
    ##geom_line(data = newdata, aes(freq_ad, fit)) +
    ##geom_ribbon(data = newdata, aes(freq_ad, ymin=lower, ymax =upper), alpha =0.2) +
    ##facet_grid(Hloc.class ~ freq.class) +
           scale_x_log10() +
    theme_bw()
## Grafico com altura em X
    ggplot(abundants3, aes(Hloc)) +
        geom_point(aes(y=ssl_tot40, color=zooc)) +
        ##geom_line(data = newdata, aes(freq_ad, fit)) +
        ##geom_ribbon(data = newdata, aes(freq_ad, ymin=lower, ymax =upper), alpha =0.2) +
        facet_grid(mass.class ~ freq.class) +
        theme_bw()


## Inspecionando efeitos com pacote effects
## Estima efeitos de freq_ad para valores fixos 
ab.ssl.tot40.selected.effects  <- predictorEffect("freq_ad",ab.ssl.tot40.full, focal.levels = 30, xlevels = 2)
plot(ab.ssl.tot40.selected.effects)




########### Limitação Temporal ###############
#O modelo agora será um glm, pois as espécies são independentes, uma vez que não há repetição nos anos
#Modelo inicial para Dredge
ab.tsl.tot12.full <- glm(cbind(nm_a_tot12, nm_p_tot12) ~ mass2 + freq2 + height2 + zooc +
                    mass2:freq2 + mass2:height2 + mass2:zooc +
                    freq2:height2 + freq2:zooc +
                    height2:zooc,
                    family=binomial,
                 data=abundants2)
                

options(na.action = "na.fail")
ab.tsl.tot12.full.d <- dredge(ab.tsl.tot12.full, beta="none")

## Modelos com deltaAIC <2
subset(ab.tsl.tot12.full.d, delta < 2) 

## modelos selecionados
ab.tsl.tot12.selected <- get.models(ab.tsl.tot12.full.d, delta < 2) 

## ICs dos efeitos de cada modelo selecionado
ab.tsl.tot12.selected.IC <- lapply(ab.tsl.tot12.selected, confint)

## Resulta numa lista de ICs
ab.tsl.tot12.selected.IC

## VIF:
ab.tsl.tot12.single <- glm(cbind(nm_a_tot12, nm_p_tot12) ~ mass2 + freq2 + height2 + zooc,
                    family=binomial,
                 data=abundants2)
vif(ab.tsl.tot12.single)





