
##########################
# Support Vector Machine #
# SVM                    #
##########################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(smotefamily)

library(e1071)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

SVM1 <- function(Train, Test) {
  
  svm <- svm(formula = Clase~.,
             kernel = "radial",
             data = Train)
  
  PredTrain <- predict(object = svm, 
                       newdata = Train, 
                       type = "response")
  
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

SVM2 <- function(Train, Test) {
  
  svm.tune <- tune(svm, Clase~., 
                   data = Train, 
                   kernel = "radial", 
                   ranges = list(cost = seq(from = 1, to = 5, by = 1), 
                                 gamma = seq(from = 0.1, to = 1, by = 0.1)))
  
  svm <- svm(formula = Clase~.,
             cost = svm.tune$best.parameters[[1]], 
             gamma = svm.tune$best.parameters[[2]],
             data = Train)
  
  PredTrain <- predict(object = svm, 
                       newdata = Train, 
                       type = "response")
  
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

SVM3 <- function(Train, Test) {
  
  TrainBalanceado <- SMOTE(X = select(Train, -Clase), 
                           target = Train$Clase, 
                           K = 5,
                           dup_size = 2)$data %>% 
    mutate(Clase = factor(class)) %>% 
    select(-class)

  svm <- svm(formula = Clase~.,
             kernel = "radial",
             data = TrainBalanceado)
    
  PredTrain <- predict(object = svm, 
                       newdata = TrainBalanceado, 
                       type = "response")
    
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
    
  MC.Train <- table(PredTrain, TrainBalanceado$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


M1.SVM <- GenerarResultadosParalelo(Metodo = "SVM1", 
                                    workers = availableCores())

M2.SVM <- GenerarResultadosParalelo(Metodo = "SVM2", 
                                    workers = availableCores())

M3.SVM <- GenerarResultadosParalelo(Metodo = "SVM3", 
                                    workers = availableCores())






