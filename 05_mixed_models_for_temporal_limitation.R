library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(lme4)
library(bbmle)
## library(rptR)
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
## Exploratory plots ##
################################################################################
## Ordered tsl
ab.sp.rsl %>%
    arrange(rank.tsl) %>%
    ggplot(aes(rank.tsl, tsl.mean)) +
    geom_point() +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max)) +
    scale_x_continuous(breaks=1:nrow(ab.sp.rsl), labels=ab.sp.rsl$species[order(ab.sp.rsl$rank.tsl)]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1) )

## Mean and range of tsl values as a function of the predictor variables    
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
## Full model to feed the dredge function
ab.tl.full <- glmer(cbind(nm_a, nm_p) ~ log.mass2 + freq2 + height2 +
                    mass2:freq2 + mass2:height2 +
                    freq2:height2 + (1|species),
                    family=binomial,
                 data=abundants,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
                                      
## Model selection with dredging
options(na.action = "na.fail")
ab.tl.full.d <- dredge(ab.tl.full, beta="none")

## Modes with deltaAIC <2
subset(ab.tl.full.d, delta<2) 

## Selected models
ab.tl.selected <- get.models(ab.tl.full.d, subset=delta<2)

############### Average model ################
tl.mavg <- model.avg(ab.tl.selected)

## ## Confidence interval of cofficients of the average model
confint(tl.mavg, full = TRUE) # full averaging 
confint(tl.mavg, full = FALSE) # conditional averaging

## Bootstrap estimation of prediction confidence intervals of the model
## conditional and unconditional to random effects.
## Data frame with data to estimate predictions
new.data.tsl <- expand.grid(
    freq2 = quantile(abundants$freq2, c(0.25, 0.75)) ,
    height2 = quantile(abundants$height2, c(0.25, 0.75)),
    log.mass2 = quantile(abundants$log.mass2, seq(0,1, by=0.05)),
    species = unique(abundants$species)[1]
)

## Function to run on each bootstrap sample
f1 <- function(.) predict(., newdata = new.data.tsl)

## Runs the bootstrap and returns a list with prediction CIs
tsl.bootMer <- ic.bootMer(lista.modelos = ab.tl.selected,
                          newdata = new.data.tsl,
                          nsim = 1000, parallel=TRUE,
                          ncpus = 6)
## Boostraped predited values
tsl.pred <- tsl.bootMer$predicted[, -4]

## Returns standardized values of predictor variables to the original scale
## (for plots)
tsl.pred %<>%
    mutate(height = height2*sd.Hloc + mean.Hloc,
           freq = freq2*sd.freq + mean.freq,
           log.mass = log.mass2*sd.log.mass + mean.log.mass,
           mass = exp(log.mass),
           pfit = ilogit(fit),
           plower = ilogit(lower),
           pupper = ilogit(upper),
           height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq)),
           )

## Plot: observed and predicted values with CI's for fixed effects and fixed + random effects
## Predicted were calculated for the median values of frequency and tree height
## of each group depicted in panels
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
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

## Logit scale
p2 <- 
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
    ggplot(aes(mass, l.tsl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin= l.tsl.upper, ymax = l.tsl.lower, color=species)) +
    geom_line(aes(y = fit), data = tsl.pred)+
    geom_ribbon(aes(y = fit, ymin = lower, ymax = upper), data = tsl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("Logit TSL")
p2

## Only the predicted lines, to evaluate effects
## Logit scale
p3 <- tsl.pred %>%
    mutate(classe = paste0(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("Predicted TSL (logit)")
p3

## Probability scale
p4 <- tsl.pred %>%
    mutate(classe = paste0(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("Predicted TSL")
p4

## All plots 
p1 ## error bars are min-max
p2 ## error bars are 1 sd ate the logit scale
p3
p4

######### Replicability (pseudo-R-squared) for the selected models ##############
## Null model to be used as a refernce for replicability calculations
ab.tl.null <- glmer(cbind(nm_a, nm_p) ~ (1|species),
                    family=binomial,
                    data=abundants)
## replicabilities, delta method
for(i in 1:length(ab.tl.selected))
    print(r.squaredGLMM(ab.tl.selected[[i]], null = ab.tl.null))
