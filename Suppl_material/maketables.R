library(knitr)
library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(lme4)
library(bbmle)
library(MuMIn)
library(modelsummary)
library(kableExtra)

## First run scripts 01 - 07 in the parent directory
## Mixed models tables as independent tex files
## Funcao para extrair coeficientes full averaged dos objetos do modelo médio
tidy_custom.averaging <- function(x, ...) {
  s <- summary(x)$coefmat.full[-1,]
  out <- data.frame(
    term = row.names(s),
    estimate = s[, "Estimate"],
    std.error = s[,"Std. Error"])
  out
}

## Funcao para incluir AICc na tabela (https://vincentarelbundock.github.io/modelsummary/articles/modelsummary.html#custom-appearance-2)
glance_custom.glmerMod <- function(x) {
    ret <- data.frame(
      AICc = AICc(x))
    ret
}
cm <- c("freq2" = "Adult Frequency", "height2" = "Adult height", "log.mass2" = "Log seed mass",
        "freq2:height2" = "Frequency:Height", "freq2:log.mass2" = "Frequency:Mass",
        "SD (Intercept)" = "SD (Intercept)", "SD (Observations)" = "SD (Observations)")
lista <- c(ab.sl.selected, list(Average = sl.mavg))
modelsummary(lista,
             gof_omit= "^(?!AICc)",
             coef_map = cm,
             output= "ssl_mixed.tex",
             title = "Mixed-effect models for the spatial seed limation. For each model is shown standardized coefficients and standard error (brackets). Also shown the estimated standard deviation for the random effects (SD), the AICc value for each selected model, and the coefficients of the average model.")

modelsummary(c(ab.tl.selected, list(Average = tl.mavg)),
             gof_omit= "^(?!AICc)",
             coef_map = cm,
             output= "tsl_mixed.tex",
             title = "Mixed-effect models for the temporal seed limation. For each model is shown standardized coefficients and standard error (brackets). Also shown the estimated standard deviation for the random effects (SD), the AICc value for each selected model, and the coefficients of the average model.")

## glms
## Function to add QAICc to the table
glance_custom.glm <- function(x) {
    ret <- data.frame(
      QAICc = QAICc(x, chat = dfun(ab.ssl.tot40.full.q)))
    ret
}

modelsummary(c(ab.ssl.tot40.selected.q, list(Average = ab.ssl.tot40.mavg)),
             gof_omit = "^(?!QAI)",
             coef_map = cm,
             output= "ssl_glm.tex",
             title = "Fixed effects models for the spatial seed limitation pooled over years (see text for details). For each model is shown standardized coefficients and standard error (brackets). Also shown  the QAICc value for each selected model, and the coefficients of the average model.")

glance_custom.glm <- function(x) {
    ret <- data.frame(
      QAICc = QAICc(x, chat = dfun(ab.tsl.tot12.full.q)))
    ret
}

modelsummary(c(ab.tsl.tot12.selected.q, list(Average = ab.tsl.tot12.mavg)),
             gof_omit = "^(?!QAI)",
             coef_map = cm,
             output= "tsl_glm.tex",
             title = "Fixed effects models for the temporal seed limitation pooled over years (see text for details). For each model is shown standardized coefficients and standard error (brackets). Also shown  the QAICc value for each selected model, and the coefficients of the average model.")
             
