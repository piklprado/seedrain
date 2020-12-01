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


################################################################################
## Mixed-effect models for spatial limitation
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

################################################################################
### Ajuste do modelo, modelo médio e tudo mais ##
################################################################################
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
           height.class = ifelse(height <= median.Hloc, paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq, paste0("Occupancy < ", median.freq), paste0("Occupancy > ", median.freq)),
           )

## O plot final: observados, previsto pelo modelo geral e IC dos fixos e total
## As especies estao divididas em 4 grupos, delimitados pelas medianas de frequencia e altura
## Os previstos são so calculados para a mediana de freq e altura em cada grupo
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste0("Height < ", median.Hloc), paste0("Height > ", median.Hloc)),
           freq.class = ifelse(freq <= median.freq, paste0("Freq < ", median.freq), paste0("Freq > ", median.freq))) %>%
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
p3 <-
    ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("SSL previsto (logito)")
p3

## Escala prob
p4 <-
    ssl.pred %>%
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
