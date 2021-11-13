library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(rmcorr)
source("funcoes.R")

#####################################################################################################################
## Total number of records in seed collectors and sampling years for all 82 species
#####################################################################################################################

## Number of seed collectors each species was recorded
par(mfrow=c(1,4), bty="n", cex.lab=1.5, mar=c(5,5,1,0))
boxplot (tudo82spp$nt_p1, ylab="Number of traps", xlab="Year 1", bty = "l", col=NULL, width=2)
boxplot (tudo82spp$nt_p2, ylab=NULL, xlab="Year 2", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nt_p3, ylab=NULL, xlab="Year 3", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nt_p_mean, ylab=NULL, xlab="Mean", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
par(mfrow=c(1,1), bty="o")

##Number of months each species was recorded
par(mfrow=c(1,4), bty="n", cex.lab=1.5, mar=c(5,5,1,0))
boxplot (tudo82spp$nm_p1, ylab="Number of months", xlab="Year 1", bty = "l", col=NULL, width=2)
boxplot (tudo82spp$nm_p2, ylab=NULL, xlab="Year 2", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nm_p3, ylab=NULL, xlab="Year 3", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
boxplot (tudo82spp$nm_p_mean, ylab=NULL, xlab="Mean", bty = "n", col=NULL, frame.plot=FALSE, yaxt="n", width=2)
par(mfrow=c(1,1), bty="o")

## Number of seed collectors over the three sampling years
boxplot (tudo82spp$nt_p_tot40, ylab="Number of traps", xlab="Total", bty = "l", col=NULL, width=2)

## Number of the 12 months of the year at which each species was recorded
boxplot (tudo82spp$nm_p_tot12, ylab="Number of months", xlab="Total", bty = "l", col=NULL, width=2)

