library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(vegan)
library(ggplot2)
library(gridExtra)
source("funcoes.R")
source("01_dataprep.R")

## Permutation according to 3 null models ##

################################################################################
## All years: shuffles values across the 36 months of sampling
################################################################################

## Null model 1: seeds of each species are distributed randomly and idenpendently acrooss months
all.null.1 <- null.model(dados = sinc.ab[,-(1:2)], nrep = 10000, FUN = "seed.rand")
## Null model 2: monthly counts of seeds of each species are shuffled across months
all.null.2 <- null.model(dados = sinc.ab[, -(1:2)], nrep = 10000, FUN = "sample")
## Null model 3: monthly counts of seed of each species are randomly moved in the time series (Rosario null model)
all.null.3 <- null.model(dados = sinc.ab[, -(1:2)], nrep = 10000, FUN = "rosario")

## Bray-Curtis distances along species in time
## Observed
all.bray <- sinc.dist(sinc.ab[,-(1:2)], method = "bray")
## From the null models
all.null.1.bray <- aaply(all.null.1, 3, sinc.dist, method = "bray")
all.null.2.bray <- aaply(all.null.2, 3, sinc.dist, method = "bray")
all.null.3.bray <- aaply(all.null.3, 3, sinc.dist, method = "bray")
## Median Bray-curtis dissimilarity
## Obs
all.dmed <- dist.med(all.bray, FUN = "mean")
## simulated
all.null.1.dmed <- aaply(all.null.1.bray, 1, dist.med, FUN = "mean")
all.null.2.dmed <- aaply(all.null.2.bray, 1, dist.med, FUN = "mean")
all.null.3.dmed <- aaply(all.null.3.bray, 1, dist.med, FUN = "mean")
## Joins medians in a single data.frame
all.null.dmed <- data.frame(Model = rep(c("Random", "Aggregated", "Autocorrelated"),each = length(all.null.1.dmed)),
                       dmed = c(all.null.1.dmed, all.null.2.dmed, all.null.3.dmed))

## P-values
sum(all.null.1.dmed >= all.dmed)/10000
sum(all.null.2.dmed >= all.dmed)/10000
sum(all.null.3.dmed >= all.dmed)/10000


################################################################################
## Each year: shuffles values across the 12 months of sampling in each year
################################################################################
## Year 1 ##
yr1.data <- sinc.ab[sinc.ab$year_study==1,-(1:2)]
yr1.data <- yr1.data[,apply(yr1.data,2,sum)>0]
## Null model 1: seeds of each species are distributed randomly and idenpendetly acrooss months
yr1.null.1 <- null.model(dados = yr1.data, nrep = 10000, FUN = "seed.rand")
## Null model 2: monthly counts of seeds of each species are shuffled across months
yr1.null.2 <- null.model(dados = yr1.data, nrep = 10000, FUN = "sample")
## Null model 3: monthly counts of seed of each species are randomly moved in the time series (Rosario null model)
yr1.null.3 <- null.model(dados = yr1.data, nrep = 10000, FUN = "rosario")
## Bray-Curtis distances along species in time
## Observed
yr1.bray <- sinc.dist(yr1.data, method = "bray")
## From the null models
yr1.null.1.bray <- aaply(yr1.null.1, 3, sinc.dist, method = "bray")
yr1.null.2.bray <- aaply(yr1.null.2, 3, sinc.dist, method = "bray")
yr1.null.3.bray <- aaply(yr1.null.3, 3, sinc.dist, method = "bray")
## Median Bray-curtis dissimilarity
## Obs
yr1.dmed <- dist.med(yr1.bray, FUN = "mean")
## simulated
yr1.null.1.dmed <- aaply(yr1.null.1.bray, 1, dist.med, FUN = "mean")
yr1.null.2.dmed <- aaply(yr1.null.2.bray, 1, dist.med, FUN = "mean")
yr1.null.3.dmed <- aaply(yr1.null.3.bray, 1, dist.med, FUN = "mean")
## Joins medians in a single data.frame
yr1.null.dmed <- data.frame(Model = rep(c("Random", "Aggregated", "Autocorrelated"),each = length(yr1.null.1.dmed)),
                       dmed = c(yr1.null.1.dmed, yr1.null.2.dmed, yr1.null.3.dmed))
## P-values
sum(yr1.null.1.dmed >= yr1.dmed)/10000
sum(yr1.null.2.dmed >= yr1.dmed)/10000
sum(yr1.null.3.dmed >= yr1.dmed)/10000

## Year 2 ##
yr2.data <- sinc.ab[sinc.ab$year_study==2,-(1:2)]
yr2.data <- yr2.data[,apply(yr2.data,2,sum)>0]
## Null model 1: seeds of each species are distributed randomly and idenpendetly acrooss months
yr2.null.1 <- null.model(dados = yr2.data, nrep = 10000, FUN = "seed.rand")
## Null model 2: monthly counts of seeds of each species are shuffled across months
yr2.null.2 <- null.model(dados = yr2.data, nrep = 10000, FUN = "sample")
## Null model 3: monthly counts of seed of each species are randomly moved in the time series (Rosario null model)
yr2.null.3 <- null.model(dados = yr2.data, nrep = 10000, FUN = "rosario")
## Bray-Curtis distances along species in time
## Observed
yr2.bray <- sinc.dist(yr2.data, method = "bray")
## From the null models
yr2.null.1.bray <- aaply(yr2.null.1, 3, sinc.dist, method = "bray")
yr2.null.2.bray <- aaply(yr2.null.2, 3, sinc.dist, method = "bray")
yr2.null.3.bray <- aaply(yr2.null.3, 3, sinc.dist, method = "bray")
## Median Bray-curtis dissimilarity
## Obs
yr2.dmed <- dist.med(yr2.bray, FUN = "mean")
## simulated
yr2.null.1.dmed <- aaply(yr2.null.1.bray, 1, dist.med, FUN = "mean")
yr2.null.2.dmed <- aaply(yr2.null.2.bray, 1, dist.med, FUN = "mean")
yr2.null.3.dmed <- aaply(yr2.null.3.bray, 1, dist.med, FUN = "mean")
## Joins medians in a single data.frame
yr2.null.dmed <- data.frame(Model = rep(c("Random", "Aggregated", "Autocorrelated"),each = length(yr2.null.1.dmed)),
                       dmed = c(yr2.null.1.dmed, yr2.null.2.dmed, yr2.null.3.dmed))
## P-values
sum(yr2.null.1.dmed >= yr2.dmed)/10000
sum(yr2.null.2.dmed >= yr2.dmed)/10000
sum(yr2.null.3.dmed >= yr2.dmed)/10000


## Year 3 ##
yr3.data <- sinc.ab[sinc.ab$year_study==3,-(1:2)]
yr3.data <- yr3.data[,apply(yr3.data,2,sum)>0]
## Null model 1: seeds of each species are distributed randomly and idenpendetly acrooss months
yr3.null.1 <- null.model(dados = yr3.data, nrep = 10000, FUN = "seed.rand")
## Null model 2: monthly counts of seeds of each species are shuffled across months
yr3.null.2 <- null.model(dados = yr3.data, nrep = 10000, FUN = "sample")
## Null model 3: monthly counts of seed of each species are randomly moved in the time series (Rosario null model)
yr3.null.3 <- null.model(dados = yr3.data, nrep = 10000, FUN = "rosario")
## Bray-Curtis distances along species in time
## Observed
yr3.bray <- sinc.dist(yr3.data, method = "bray")
## From the null models
yr3.null.1.bray <- aaply(yr3.null.1, 3, sinc.dist, method = "bray")
yr3.null.2.bray <- aaply(yr3.null.2, 3, sinc.dist, method = "bray")
yr3.null.3.bray <- aaply(yr3.null.3, 3, sinc.dist, method = "bray")
## Median Bray-curtis dissimilarity
## Obs
yr3.dmed <- dist.med(yr3.bray, FUN = "mean")
## simulated
yr3.null.1.dmed <- aaply(yr3.null.1.bray, 1, dist.med, FUN = "mean")
yr3.null.2.dmed <- aaply(yr3.null.2.bray, 1, dist.med, FUN = "mean")
yr3.null.3.dmed <- aaply(yr3.null.3.bray, 1, dist.med, FUN = "mean")
## Joins medians in a single data.frame
yr3.null.dmed <- data.frame(Model = rep(c("Random", "Aggregated", "Autocorrelated"),each = length(yr3.null.1.dmed)),
                       dmed = c(yr3.null.1.dmed, yr3.null.2.dmed, yr3.null.3.dmed))
## P-values
sum(yr3.null.1.dmed >= yr3.dmed)/10000
sum(yr3.null.2.dmed >= yr3.dmed)/10000
sum(yr3.null.3.dmed >= yr3.dmed)/10000

################################################################################
## Plots
################################################################################

## Histograms of Bray-curtis dissimilarities for a quick inspection
hist(all.bray[lower.tri(all.bray)])
abline(v = dist.med(all.bray, FUN= "median"))
abline(v = dist.med(all.bray, FUN= "mean"), col ="blue")

## hist(yr1.bray)
## hist(yr2.bray)
## hist(y3.bray)

## Density plots of null distributions and observed values
p <- all.null.dmed %>%
    ggplot(aes(dmed, group = Model)) +
    geom_density(aes(color = Model), alpha = 0.5) +
    scale_fill_manual(values = c("green", "red", "blue")) +
    xlab("Median Temporal Dissimilarity (Bray-Curtis)") +
    theme_classic()

p.all <-
    p +
    geom_vline(xintercept = all.dmed, lty =2, col = "darkblue") +
    ggtitle("All Years") +
    theme(legend.position = c(0.2, 0.75))

p1 <- p %+% yr1.null.dmed  +
    geom_vline(xintercept = yr1.dmed, lty =2, col = "darkblue") +
    theme(legend.position = "none") +
    ggtitle("Year 1")

p2 <- p %+% yr2.null.dmed  +
    geom_vline(xintercept = yr2.dmed, lty =2, col = "darkblue") +
    theme(legend.position = "none") +
    ggtitle("Year 2")

p3 <- p %+% yr3.null.dmed  +
    geom_vline(xintercept = yr3.dmed, lty =2, col = "darkblue") +
    theme(legend.position = "none") +
    ggtitle("Year 3")

grid.arrange(p.all, p1, p2, p3,
             nrow = 2, ncol = 2)

## Plot of points & errorbars of all null models
## Preparing dataframes
summary.dmed <- cbind(
    year = rep(c("All Years", "Year 1", "Year 2", "Year 3"), each = 3),
    rbind(
        group_by(all.null.dmed, by = Model) %>% summarise(mean = mean(dmed),
                                                          mean.low = quantile(dmed, 0.025),
                                                          mean.upp = quantile(dmed, 0.975)) %>% as.data.frame(),
        group_by(yr1.null.dmed, by = Model) %>% summarise(mean = mean(dmed),
                                                          mean.low = quantile(dmed, 0.025),
                                                          mean.upp = quantile(dmed, 0.975)) %>% as.data.frame(),
        group_by(yr2.null.dmed, by = Model) %>% summarise(mean = mean(dmed),
                                                          mean.low = quantile(dmed, 0.025),
                                                          mean.upp = quantile(dmed, 0.975)) %>% as.data.frame(),
        group_by(yr3.null.dmed, by = Model) %>% summarise(mean = mean(dmed),
                                                          mean.low = quantile(dmed, 0.025),
                                                          mean.upp = quantile(dmed, 0.975)) %>% as.data.frame()
        )
)
summary.dmed$by <- factor(summary.dmed$by, levels=c("Random", "Aggregated", "Autocorrelated"))
summary.dmed %<>% mutate(mean.ov = mean, mean.ov.low = mean.low, mean.ov.upp = mean.upp)

observed <- data.frame(year = c("All Years", "Year 1", "Year 2", "Year 3"),
                       by = "Observed", mean = c(all.dmed, yr1.dmed, yr2.dmed, yr3.dmed))
observed %<>% mutate(mean.ov = mean)

## The plot
## Colorblind-friendly pallete (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png("figures/null_models_synchr.png")
ggplot(summary.dmed, aes(year, mean.ov)) +
    geom_point(data = observed, aes(y = mean.ov, shape = "Observed"),
               size = 3, position=position_nudge(x = -0.3)) +
    geom_errorbar(aes(ymin = mean.ov.low, ymax = mean.ov.upp, color = by),
                  position = position_dodge(width=0.2), size = 1.5, alpha =0.75) +
    ylab("Mean temporal asynchrony (Bray Curtis)") +
    xlab("") +
    labs(color = "Null model:", shape = "") +
    scale_color_manual(values = cbPalette[c(2,3,7)]) +
    theme_classic() +
    theme(axis.text = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1.25)),
          legend.position = "bottom")
dev.off()    
    
