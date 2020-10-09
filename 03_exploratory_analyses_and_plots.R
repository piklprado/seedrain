library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(rmcorr)
source("funcoes.R")

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

