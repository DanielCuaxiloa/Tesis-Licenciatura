
#################################
# Network Discriminant Analysis #
# - Network LDA                 #
# - Network QDA                 #
#################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(ggplot2)

library(readr)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

source("../../Network Discriminant Analysis/NetQDA.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

Network.QDA <- function(Train, Test) {
  
  NetQDA <- NetQDA(formula = Clase~., 
                   datos = Train,
                   rho = 0.05)
  
  PredTrain <- Predict.NetQDA(object = NetQDA,
                              NewData = select(Train, -Clase))
  
  PredTest <- Predict.NetQDA(object = NetQDA,
                             NewData = select(Test, -Clase))
  
  MC.Train <- table(Train$Clase, PredTrain$clase_predicha$.)
  MC.Test <- table(Test$Clase, PredTest$clase_predicha$.)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


