
################################
# Linear Discriminant Analysis #
# LDA                          #
################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(smotefamily)

library(MASS)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

LDA1 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~., 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class

  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LDA2 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~.^2, 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class
  
  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LDA3 <- function(Train, Test) {
  
  TrainBalanceado <- SMOTE(X = select(Train, -Clase), 
                           target = Train$Clase, 
                           K = 5,
                           dup_size = 2)$data %>% 
    mutate(Clase = factor(class)) %>% 
    select(-class)

  LDA <- lda(formula = Clase~., 
             data = TrainBalanceado)
  
  PredTrain <- predict(object = LDA,
                       newdata = TrainBalanceado)$class
  
  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, TrainBalanceado$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


M1.LDA <- GenerarResultadosParalelo(Metodo = "LDA1", 
                                    workers = availableCores())

M2.LDA <- GenerarResultadosParalelo(Metodo = "LDA2", 
                                    workers = availableCores())

M3.LDA <- GenerarResultadosParalelo(Metodo = "LDA3", 
                                    workers = availableCores())


