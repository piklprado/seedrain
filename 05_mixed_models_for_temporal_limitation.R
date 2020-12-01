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
source("01_dataprep.R")

################################################################################
## Mixed-Effects models for temporal limitation
################################################################################

################################################################################
## Graficos exploratorios ##
################################################################################
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

################################################################################
## Model fit and model averaging
################################################################################
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


## Calculo do intervalo de previsao com efeito fixo e aleatorio, por
## bootstrap
## Primeiro criamos uma dataframe com os valores de altura a
## prever e uma espécie qualquer 
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
           height.class = ifelse(height <= median.Hloc, paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq, paste0("Occupancy < ", median.freq), paste0("Occupancy > ", median.freq)),
           )

## O plot final: observados, previsto pelo modelo geral e IC dos fixos e total
## As especies estao divididas em 4 grupos, delimitados pelas medianas de frequencia e altura
## Os previstos são so calculados para a mediana de freq e altura em cada grupo
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq, paste0("Occupancy < ", median.freq), paste0("Occupancy > ", median.freq))) %>%
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
    mutate(height.class = ifelse(height <= median.Hloc, paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq, paste0("Occupancy < ", median.freq), paste0("Occupancy > ", median.freq))) %>%
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
    mutate(classe = paste0(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("TSL previsto (logito)")
p3

## Escala prob
p4 <- tsl.pred %>%
    mutate(classe = paste0(height.class,freq.class, sep =" , ")) %>%
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
