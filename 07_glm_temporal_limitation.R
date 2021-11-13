library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(MuMIn)
library(car)
source("funcoes.R")
source("01_dataprep.R")

################################################################################
##  Model fit
################################################################################
## Full model to feed dredge function (zoocoric species excluded)
## Binomial fit
ab.tsl.tot12.full <- glm(cbind(nm_a_tot12, nm_p_tot12) ~ log.mass2 + freq2 + height2 +
                             log.mass2:freq2 + log.mass2:height2 + freq2:height2,
                         family="binomial",
                         data=abundants2)
## Quasibinomial fit
## Tweak to allow dredge work with quasibinomial
ab.tsl.tot12.full.q <- update(ab.tsl.tot12.full, family="x.quasibinom")

## Dredge
options(na.action = "na.fail")
## Binomial
ab.tsl.tot12.full.d <- dredge(ab.tsl.tot12.full, beta="none", rank=AICc)
## Quasibinomial
ab.tsl.tot12.full.q.d <- dredge(ab.tsl.tot12.full.q, beta="none", rank = "QAICc", chat = dfun(ab.tsl.tot12.full.q))

## Models with deltaAIC <2 ##
# Binomial
subset(ab.tsl.tot12.full.d, delta < 2)
## quasibinomial
subset(ab.tsl.tot12.full.q.d, delta < 2) 

## Selected models
## Binomial
ab.tsl.tot12.selected <- get.models(ab.tsl.tot12.full.d, delta < 2)
## Quasibinomial
ab.tsl.tot12.selected.q <- get.models(ab.tsl.tot12.full.q.d, delta < 2)

## CIs of the coefficients of each selected model
## Binomial
ab.tsl.tot12.selected.IC <- lapply(ab.tsl.tot12.selected, confint)
## Quasibinomial
ab.tsl.tot12.selected.IC.q <- lapply(ab.tsl.tot12.selected.q, confint)

## VIF: 
ab.tsl.tot12.single <- glm(cbind(nm_a_tot12, nm_p_tot12) ~ log.mass2 + freq2 + height2 + zooc,
                    family=binomial, data=abundants2)
vif(ab.tsl.tot12.single)

## Pseudo R-squared
ab.tsl.tot12.null <- glm(cbind(nm_a_tot12, nm_p_tot12) ~ 1,
                         family="x.quasibinom",
                         data=abundants2)
for(i in 1:length(ab.tsl.tot12.selected.q))
    print(r.squaredGLMM(ab.tsl.tot12.selected.q[[i]], null = ab.tsl.tot12.null))

################################################################################
## Model averaging
################################################################################
ab.tsl.tot12.mavg <- model.avg(ab.tsl.tot12.selected.q)

## CIs of coefficients of average model
confint(ab.tsl.tot12.mavg, full=TRUE) # full
confint(ab.tsl.tot12.mavg, full=FALSE) # conditional

## Predicted values 
## Dataframe to make the predictions
ab.tsl.tot12.newdata <- expand.grid(log.mass2 = seq(min(abundants2$log.mass2), max(abundants2$log.mass2), length =30),
                       height2 = quantile(abundants2$height2, c(0.25, 0.75)),
                       freq2 = quantile(abundants2$freq2, c(0.25, 0.75)))
## Predicted at response scale
predict <- as.data.frame(predict(ab.tsl.tot12.mavg,
                                 newdata = ab.tsl.tot12.newdata,
                                 se.fit=TRUE))
ab.tsl.tot12.newdata <- cbind(ab.tsl.tot12.newdata, predict)
## add unscaled variables to predicted values
ab.tsl.tot12.newdata$mass  <- with(abundants2, exp( (ab.tsl.tot12.newdata$log.mass2 * sd(log(mass)) + mean(log(mass)) )))
## Unscaled Cutoffs for height and frequency
Hloc.median <- median(abundants2$Hloc)
log.mass.median  <-  median(log(abundants2$mass))
freq.median  <-  median(abundants2$freq_ad)
## Adds CIs and other auxiliary variables to the dataframe 
ab.tsl.tot12.newdata  <- mutate(ab.tsl.tot12.newdata,
                                Hloc.class = ifelse(height2 < median(height2),
                                                    paste0("Tree Height < ", Hloc.median, " m"),
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
    mutate(Hloc.class = ifelse(Hloc < Hloc.median,
                               paste0("Tree Height < ", Hloc.median, " m"),
                               paste0("Tree Height > ", Hloc.median, " m")),
           mass.class = ifelse(log(mass) < log.mass.median, paste0("Seed Mass < ", exp(log.mass.median), " g"),
                               paste0("Seed Mass > ", exp(log.mass.median), " g")),
           freq.class =ifelse(freq_ad < freq.median, paste0("Occupancy < ", freq.median), paste0("Occupancy > ", freq.median)) ,
           ssl_tot12 = nm_a_tot12/12, 
           logito.ssl = log(nm_a_tot12/nm_p_tot12),
           tsl_tot12 = nm_a_tot12/12, 
           logito.tsl = log(nm_a_tot12/nm_p_tot12))
## Observed x predicted plots ##
## Plot in response scale
ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=tsl_tot12, color=zooc)) +
    geom_line(data = ab.tsl.tot12.newdata, aes(mass, fit.resp)) +
    geom_ribbon(data = ab.tsl.tot12.newdata, aes(mass, ymin=lower.resp, ymax =upper.resp), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
    scale_x_log10() +
    theme_bw()
## Plot in logit scale
ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=logito.tsl, color=zooc)) +
    geom_line(data = ab.tsl.tot12.newdata, aes(mass, fit)) +
    geom_ribbon(data = ab.tsl.tot12.newdata, aes(mass, ymin=lower, ymax =upper), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
           scale_x_log10() +
    theme_bw()


