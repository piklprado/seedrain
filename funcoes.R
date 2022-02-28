#' logit function
#'
logito <- function(x, max =1, delta =0){
    log(x/(max-x))
}

#' Inverse logit 
#' @param x logit 
ilogit <- function(x)
    exp(x)/(1+exp(x))

#' Prediction interval for average model from glmms
#' @param lista.modelos lis of model objects that make the average model
#' @param newdata new data to make the predictions
#' @param nsim number of simulations
#' @param ... further arguments to the bootMer function
ic.bootMer <- function(lista.modelos, newdata, nsim = 1000, ...){
    pesos <- AICctab(lista.modelos, weights = TRUE, mnames = names(lista.modelos))$weight
    nsims <- round(nsim * pesos)
    f1 <- function(.) predict(., newdata = newdata)
    funcao.boot <- function(x, N)
       bootMer(x, FUN = f1, nsim = N, ...)$t 
    lista.boot <- mapply(funcao.boot, x = lista.modelos, N = nsims)
    df.boot <- ldply(lista.boot)
    newdata$fit  <- apply(df.boot[,-1], 2, mean)
    newdata$lower <- apply(df.boot[,-1], 2, quantile, 0.025)
    newdata$upper <- apply(df.boot[,-1], 2, quantile, 0.975)
    list(predicted = newdata, boot = lista.boot, boot.t = df.boot)
}

#' Randomize seeds over months
#'
seed.rand <- function(x){
    indexes <- sample(factor(1:length(x)), size = sum(x), replace = TRUE)
    as.numeric(table(indexes))
}

#' Shuffles times series using Rosario null model
#' @param x a vector of values in temporal sequence
#' @references Castro‐Arellano, I., Lacher Jr, T. E., Willig, M. R., &
#'     Rangel, T. F. (2010). Assessment of assemblage‐wide temporal
#'     niche segregation using null models. Methods in ecology and
#'     evolution, 1(3), 311-318.
rosario <- function(x){
    L <- length(x)
    indices <- 1:L
    shift <- sample(0:L, size=1)
    rand.ind <- (indices + shift)%%L
    rand.ind[rand.ind==0] <- L
    if(runif(1)>0.5)
        rand.ind <- rev(rand.ind)
    return(x[rand.ind])
    }

#' Shuffles a matrix of month x species and outuputs similarities
#' matrices
#' @param dados: matrix or dataframe with abundances of species
#'     (columns) by mnonths (rows)
#' @param estrato vector with the variable to stratify the
#'     randomization. The dataset will be subdivided by the levels of
#'     this factor, and then randomization is applied to each of these
#'     subdivisions. If absent, randomization will be applied to the
#'     dataset without splitting it.
#' @param nrep number of repetitions of the randomization
#' @param FUN null model to be used. "seed.rand" distributes seeds
#'     individually and independently of each species between
#'     months. "sample" keeps the number of seeds in each month of
#'     each species, but shuffles the months. "rosario" keeps the
#'     total number of seeds per month and the succession of values
#'     ​​over the time series, only changing the position of the series
#'     over time (Castro-Arellano et al 2010)
#' @return an array with dimension [nrow(data), ncol(data), nrep];
#'     with data matrices whose columns were permuted by the null
#'     model, with or without stratification.
#' @references Castro‐Arellano, I., Lacher Jr, T. E., Willig, M. R., &
#'     Rangel, T. F. (2010). Assessment of assemblage‐wide temporal
#'     niche segregation using null models. Methods in ecology and
#'     evolution, 1(3), 311-318.
null.model <- function(dados, estrato, nrep = 100, FUN = c("seed.rand", "sample", "rosario")){
    f2 <- match.fun(FUN)
    if(missing(estrato))
        f1 <- function(x) as.matrix(apply(x, 2, f2))
    else
        f1 <- function(x) as.matrix(ddply(x, .(estrato), .fun = function(y) apply(y, 2, f2))[,-1])
    results <- array(dim=c(nrow(dados), ncol(dados), nrep))
    for(i in 1:nrep)
        results[,,i] <- f1(dados)
    results
}



#' Distances to the n-th nearest neigbor
#' @param x distance matrix
#' @return a data frame with the n-th distance ranking of each neighbor and the respective distance.
dvmp <- function(x){
    m1 <- apply(x, 2, sort)[-1,]
    data.frame(ranking = 1:nrow(m1),
               media = apply(m1, 1, mean),
               mediana = apply(m1, 1, median)
               ##lower = apply(m1, 1, quantile, 0.025),
               ##upper = apply(m1, 1, quantile, 0.975)
               )
    }

#'Auxiliary functions for quasibinomial model selection and model averaging
#' @source https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf
dfun <- function(object) {
    with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}
x.quasibinom <- function(...) {
    res <- quasibinomial(...)
    res$aic <- binomial(...)$aic
    res
}

#' Auxiliary function to return a distance matrix as a matrix
sinc.dist <- function(x, ...)
    as.matrix(vegdist(t(x), ...))

#' Auxiliary function to calculate a summary function from distance matrices
dist.med <- function(x, FUN = "mean"){
    f1 <- match.fun(FUN)
    f1(x[lower.tri(x)])
}

#' Auxiliary function: relative abundance
ab.rel <- function(x){
    if(all(x==0))
        x
    else
        x/sum(x)
}
