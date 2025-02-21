
#################################
# Recopilación de Resultados TE #
#################################

library(dplyr)
library(forcats)
library(ggplot2)

Modelo1 <- read.csv(file = "TE_Modelo1.csv", 
                    stringsAsFactors = TRUE)

Modelo2 <- read.csv(file = "TE_Modelo2.csv", 
                    stringsAsFactors = TRUE)

Modelo3 <- read.csv(file = "TE_Modelo3.csv", 
                    stringsAsFactors = TRUE)

Modelo4 <- read.csv(file = "TE_Modelo4.csv", 
                    stringsAsFactors = TRUE)

Modelo5 <- read.csv(file = "TE_Modelo5.csv", 
                    stringsAsFactors = TRUE)

Datos <- bind_rows(Modelo1, Modelo2, Modelo3, Modelo4, Modelo5)
