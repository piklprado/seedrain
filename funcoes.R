#' funcao logito
#'
logito <- function(x, max =1, delta =0){
    log(x/(max-x))
}

#' Inverso da funcao logito
#' @param x logito 
ilogit <- function(x)
    exp(x)/(1+exp(x))

#' Estima intervalos de previsao de modelo medio com efeitos mistos
#' @param lista.modelos lista de modelos dos quais se deseja fazer o
#'     modelo medio
#' @param newdata planilha com os valores das preditoras para os quais
#'     se deseja fazer as previsoes (veja argumento newdata em
#'     predict)
#' @param nsim total de simulacoes a fazer
#' @param ... outros argumentos a passar para a função bootMer
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

#' Aleatorização de series temporais pelo algoritmo Rosario
#' @details aplica a um vetor a permutação circular "rosario", criada
#'     para séries temporais
#' @param x um vetor
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

#' Aleatoriza a matriz de meses x especies e devolve as matrizes de
#' similaridade
#' @param dados: matriz ou dataframe apenas com as abundancias das
#'     especies (colunas) por meses(linhas)
#' @param estrato vetor com a variavel para estratificar a
#'     randomização. O conjunto de dados será subdividido pelos níveis
#'     destes estrato, e então a aleatorização aplicada a cada uma
#'     destas subdivisões. Se ausente, a aleatorização será aplicada
#'     ao conjunto de dados sem dividí-lo.
#' @param nrep número de repetições da randomização
#' @param FUN algoritmo de modelos nulo a utilizar. "seed,rand"
#'     distribui sementes individualmente e independentemente de cada
#'     espécies entre os meses. "sample" mantém o número de sementes
#'     em cada mês de cada espécie, mas embaralha os meses. "rosario"
#'     mantém o total o número de sementes por mês e a sucessão de
#'     valores ao longo da série temporal, apenas alterando a posição
#'     da série ao longo do tempo (Castro-Arellano et al 2010)
#' @return um array de dimensão [nrow(dados), ncol(dados), nrep]; com
#'     matrizes de dados cujas colunas forma permutadas pelo algoritmo
#'     rosario, com ou sem estratificação.
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



#' Distancias ao n-esimo vizinho mais próximo de uma matriz de distâncias
#' @param x matriz de distâncias
#' @return data.frame com ranking do vizinho, distância média e limites do IC a 95%
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
## Funcao para aplicar a cada ano: abundancia relativa
ab.rel <- function(x){
    if(all(x==0))
        x
    else
        x/sum(x)
}
