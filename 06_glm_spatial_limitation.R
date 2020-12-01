library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ciTools)
library(ggplot2)
library(MuMIn)
library(AICcmodavg)
library(effects)
library(car)
source("funcoes.R")
source("01_dataprep.R")

################################################################################
##  Model fit
################################################################################
## Modelo inicial para Dredge (SEM ZOOCORIA)
## Binomial fit
ab.ssl.tot40.full <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ log.mass2 + freq2 + height2 +
                             log.mass2:freq2 + log.mass2:height2 + freq2:height2,
                         family="binomial",
                         data=abundants2)
## Quasibinomial fit
## Modelo ajustado de tal forma que o dredge funfe com quasibinomial
ab.ssl.tot40.full.q <- update(ab.ssl.tot40.full, family="x.quasibinom")

## Dredge
options(na.action = "na.fail")
## Binomial
ab.ssl.tot40.full.d <- dredge(ab.ssl.tot40.full, beta="none", rank=AICc)
## Quasibinomial
ab.ssl.tot40.full.q.d <- dredge(ab.ssl.tot40.full.q, beta="none", rank = "QAICc", chat = dfun(ab.ssl.tot40.full.q))

## Modelos com deltaAIC <2 ##
# Binomial
subset(ab.ssl.tot40.full.d, delta < 2)
## quasibinomial
subset(ab.ssl.tot40.full.q.d, delta < 2) 

## modelos selecionados
## Binomial
ab.ssl.tot40.selected <- get.models(ab.ssl.tot40.full.d, delta < 2)
## Quasibinomial
ab.ssl.tot40.selected.q <- get.models(ab.ssl.tot40.full.q.d, delta < 2)

## ICs dos efeitos de cada modelo selecionado
## Binomial
ab.ssl.tot40.selected.IC <- lapply(ab.ssl.tot40.selected, confint)
## Quasibinomial
ab.ssl.tot40.selected.IC.q <- lapply(ab.ssl.tot40.selected.q, confint)

## VIF: 
ab.ssl.tot40.single <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ log.mass2 + freq2 + height2 + zooc,
                    family=binomial, data=abundants2)

vif(ab.ssl.tot40.single)

## Pseudo-Rsquared
ab.ssl.tot40.null <- glm(cbind(nt_a_tot40, nt_p_tot40) ~ 1,
                         family="x.quasibinom",
                         data=abundants2)
for(i in 1:length(ab.ssl.tot40.selected.q))
    print(r.squaredGLMM(ab.ssl.tot40.selected.q[[i]], null = ab.ssl.tot40.null))

################################################################################
## Model averaging
################################################################################
ab.ssl.tot40.mavg <- model.avg(ab.ssl.tot40.selected.q)

## CIs dos efeitos medios
confint(ab.ssl.tot40.mavg, full=TRUE) # para full
confint(ab.ssl.tot40.mavg, full=FALSE) # para full

## Previstos
## Dataframe to make the predictions
ab.ssl.tot40.newdata <- expand.grid(log.mass2 = seq(min(abundants2$log.mass2), max(abundants2$log.mass2), length =30),
                       height2 = quantile(abundants2$height2, c(0.25, 0.75)),
                       freq2 = quantile(abundants2$freq2, c(0.25, 0.75)))
## Predicted at response scale
predict <- as.data.frame(predict(ab.ssl.tot40.mavg,
                                 newdata = ab.ssl.tot40.newdata,
                                 se.fit=TRUE))
ab.ssl.tot40.newdata <- cbind(ab.ssl.tot40.newdata, predict)
## add unscaled variables to predicted values
ab.ssl.tot40.newdata$mass  <- with(abundants2, exp( (ab.ssl.tot40.newdata$log.mass2 * sd(log(mass)) + mean(log(mass)) )))
## Unscaled Cutoffs for height and frequency
Hloc.median <- median(abundants2$Hloc)
log.mass.median  <-  median(log(abundants2$mass))
freq.median  <-  median(abundants2$freq_ad)
## Adds CIs and other auxiliary variables to the dataframe 
ab.ssl.tot40.newdata  <- mutate(ab.ssl.tot40.newdata,
                                Hloc.class = ifelse(height2 < median(height2), paste0("Tree Height < ", Hloc.median, " m"),
                                                    paste0("Tree Height > ", Hloc.median, " m")),
                   freq.class = ifelse(freq2 < median(freq2),
                                       paste0("Occupancy < ", freq.median), paste0("Occupancy > ", freq.median)),
                   lower = fit - 2*se.fit,
                   upper = fit + 2*se.fit,
                   fit.resp = ilogit(fit),
                   lower.resp = ilogit(lower),
                   upper.resp = ilogit(upper))
## A new dataframe for the predictions
abundants3 <- abundants2 %>%
    mutate(Hloc.class = ifelse(Hloc < Hloc.median, paste0("Tree Height < ", Hloc.median, " m"),
                               paste0("Tree Height > ", Hloc.median, " m")),
           mass.class = ifelse(log(mass) < log.mass.median, paste0("Seed Mass < ", exp(log.mass.median), " g"),
                               paste0("Seed Mass > ", exp(log.mass.median), " g")),
           freq.class =ifelse(freq_ad < freq.median, paste0("Occupancy < ", freq.median), paste0("Occupancy > ", freq.median)) ,
           ssl_tot40 = nt_a_tot40/40, 
           logito.ssl = log(nt_a_tot40/nt_p_tot40),
           tsl_tot12 = nm_a_tot12/12, 
           logito.tsl = log(nm_a_tot12/nm_p_tot12))
## Observed x predicted plots ##
## Plot in response scale
ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=ssl_tot40, color=zooc)) +
    geom_line(data = ab.ssl.tot40.newdata, aes(mass, fit.resp)) +
    geom_ribbon(data = ab.ssl.tot40.newdata, aes(mass, ymin=lower.resp, ymax =upper.resp), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
    scale_x_log10() +
    theme_bw()
## Plot in logit scale
ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=logito.ssl, color=zooc)) +
    geom_line(data = ab.ssl.tot40.newdata, aes(mass, fit)) +
    geom_ribbon(data = ab.ssl.tot40.newdata, aes(mass, ymin=lower, ymax =upper), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
           scale_x_log10() +
    theme_bw()


