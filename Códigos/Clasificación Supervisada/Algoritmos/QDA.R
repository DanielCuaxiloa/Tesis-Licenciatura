
###################################
# Quadratic Discriminant Analysis #
# QDA                             #
###################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(smotefamily)

library(MASS)


# Funciones Auxiliares ----------------------------------------------------

source("FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("Folds.RData")


# Esquemas de clasificación -----------------------------------------------

QDA1 <- function(Train, Test) {
  
  QDA <- qda(formula = Clase~., 
             data = Train)
  
  PredTrain <- predict(object = QDA,
                       newdata = Train)$class
  
  PredTest <- predict(object = QDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

QDA2 <- function(Train, Test) {
  
  QDA <- qda(formula = Clase~.^2, 
             data = Train)
  
  PredTrain <- predict(object = QDA,
                       newdata = Train)$class
  
  PredTest <- predict(object = QDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

QDA3 <- function(Train, Test) {
  
  TrainBalanceado <- SMOTE(X = select(Train, -Clase), 
                           target = Train$Clase, 
                           K = 5,
                           dup_size = 2)$data %>% 
    mutate(Clase = factor(class)) %>% 
    select(-class)
  
  QDA <- qda(formula = Clase~., 
             data = TrainBalanceado)
  
  PredTrain <- predict(object = QDA,
                       newdata = TrainBalanceado)$class
  
  PredTest <- predict(object = QDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, TrainBalanceado$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


M1.QDA <- GenerarResultadosParalelo(Metodo = "QDA1", 
                                    workers = availableCores())

M2.QDA <- GenerarResultadosParalelo(Metodo = "QDA2", 
                                    workers = availableCores())

M3.QDA <- GenerarResultadosParalelo(Metodo = "QDA3", 
                                    workers = availableCores())
