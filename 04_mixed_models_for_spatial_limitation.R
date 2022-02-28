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
## Mixed-effect models for spatial limitation
################################################################################
## Exploratory plots ##
## Ordered ssl
ab.sp.rsl %>%
    arrange(rank.ssl) %>%
    ggplot(aes(rank.ssl, ssl.mean)) +
    geom_point() +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max)) +
    scale_x_continuous(breaks=1:nrow(ab.sp.rsl), labels= ab.sp.rsl$species[order(ab.sp.rsl$rank.ssl)]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1))

## Mean and range of ssl values as a function of the predictor variables    
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
### Model fitting and model averaging ##
################################################################################
## Full model to feed the dredge function
ab.sl.full <- glmer(cbind(nt_a, nt_p) ~ log.mass2 + freq2 + height2 +
                    log.mass2:freq2 + log.mass2:height2 +
                    freq2:height2 + (1|species),
                    family=binomial,
                 data=abundants,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

## Model selection with dredging
options(na.action = "na.fail")
ab.sl.full.d <- dredge(ab.sl.full, beta="none")

## Models with  deltaAIC <2
subset(ab.sl.full.d, delta < 2) 

## Selected models
ab.sl.selected <- get.models(ab.sl.full.d, delta < 2)

## Confidence intervals of the coefficients of all selected models
ab.sl.selected.IC <- lapply(ab.sl.selected, confint)

## Checking variation inflation factors
vif(ab.sl.selected[[3]])

############### Average model ################
sl.mavg <- model.avg(ab.sl.selected)

## Confidence interval of cofficients of the average model
confint(sl.mavg, full=TRUE) # full averaging
confint(sl.mavg, full=FALSE) # conditional avergaing

## Bootstrap estimation of prediction confidence intervals of the model
## conditional and unconditional to random effects.
## Data frame with data to estimate predictions
new.data.ssl <- expand.grid(
    freq2 = quantile(abundants$freq2, c(0.25, 0.75)) ,
    height2 = quantile(abundants$height2, c(0.25, 0.75)),
    log.mass2 = quantile(abundants$log.mass2, seq(0,1, by=0.05)),
    species = unique(abundants$species)[1]
)

## Function to run on each bootstrap sample
f1 <- function(.) predict(., newdata = new.data.ssl)

## Runs the bootstrap and returns a list with prediction CIs
ssl.bootMer <- ic.bootMer(lista.modelos = ab.sl.selected,
                          newdata = new.data.ssl,
                          nsim = 1000, parallel=TRUE,
                          ncpus = 6)
## Boostraped predited values
ssl.pred <- ssl.bootMer$predicted[, -4]
## Returns standardized values of predictor variables to the original scale
## (for plots)
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
           freq.class = ifelse(freq <= median.freq, paste0("Adult Frequency < ", median.freq), paste0("Adult Frequency > ", median.freq)),
           )

## Plot: observed and predicted values with CI's for fixed effects and fixed + random effects
## Predicted were calculated for the median values of frequency and tree height
## of each group depicted in panels
p1 <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc, paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq, paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
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
    
## Same plot in logit scale
p2 <- 
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
    ggplot(aes(mass, l.ssl.mean)) +
    geom_point(aes(color=species), size=2) +
    geom_linerange(aes(ymin= l.ssl.upper, ymax = l.ssl.lower, color=species)) +
    geom_line(aes(y = fit), data = ssl.pred)+
    geom_ribbon(aes(y = fit, ymin = lower, ymax = upper), data = ssl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("Logit SSL")
p2

## Only the predicted lines, to evaluate effects
## Logit scale
p3 <-
    ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=fit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = lower, ymax =upper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw()+
    ylab("Predicted SSL (logit)")
p3

## Probability scale
p4 <-
    ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("Predicted SSL")
p4
    
## All plots
p1 ## error bars are min-max
p2 ## error bars are 1 sd ate the logit scale
p3
p4

######### Replicability (pseudo-R-squared) for the selected models ##############
## Null model to be used as a refernce for replicability calculations
ab.sl.null <- glmer(cbind(nt_a, nt_p) ~ (1|species),
                    family=binomial,
                    data=abundants)
## replicabilities, with delta method
for(i in 1:length(ab.sl.selected))
    print(r.squaredGLMM(ab.sl.selected[[i]], null = ab.sl.null))
