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

load("../.RData")
gm <- modelsummary::gof_map %>% filter(raw == "aic" )
gm$omit <- FALSE
cm <- c("freq2" = "Adult Frequency", "height2" = "Adult height", "log.mass2" = "Log seed mass",
        "freq2:height2" = "Frequency:Height", "freq2:log.mass2" = "Frequency:Mass",
        "SD (Intercept)" = "SD (Intercept)", "SD (Observations)" = "SD (Observations)")
lista <- c(ab.sl.selected, list(Average = sl.mavg))
modelsummary(lista,
             gof_map = gm,
             coef_map = cm,
             output= "ssl_mixed.tex",
             title = "Mixed-effect models for the spatial seed limation. For each model is shown standarzide coefficients and standard error (brackets).")
             
