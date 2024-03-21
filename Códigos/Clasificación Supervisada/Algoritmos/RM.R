
################################
# Multinomial Regression lasso #
# MR lasso                     #
################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(smotefamily)

library(glmnet)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

MR_lasso1 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~., 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~., 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 10,
                         type.measure = "class",
                         family = "multinomial", 
                         type.multinomial = "ungrouped")
  
  PredTrain <- predict(object = lasso.tun, 
                       newx = XTrain, 
                       type = "class", 
                       s = "lambda.min")
  
  PredTest <- predict(object = lasso.tun, 
                      newx = XTest, 
                      type = "class", 
                      s = "lambda.min")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

MR_lasso2 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 10,
                         lambda = seq(from = 0, to = 15, by = 0.1),
                         type.measure = "class",
                         family = "multinomial", 
                         type.multinomial = "ungrouped")
  
  PredTrain <- predict(object = lasso.tun, 
                       newx = XTrain, 
                       type = "class", 
                       s = "lambda.min")
  
  PredTest <- predict(object = lasso.tun, 
                      newx = XTest, 
                      type = "class", 
                      s = "lambda.min")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

MR_lasso3 <- function(Train, Test) {
  
  TrainBalanceado <- SMOTE(X = select(Train, -Clase), 
                           target = Train$Clase, 
                           K = 5,
                           dup_size = 2)$data %>% 
    mutate(Clase = factor(class)) %>% 
    select(-class)
  
  XTrainBalanceado <- model.matrix(Clase~., 
                                   data = TrainBalanceado)[,-1]
  YTrain <- TrainBalanceado$Clase
  
  XTest <- model.matrix(Clase~., 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrainBalanceado, 
                         y = YTrain, 
                         nfolds = 10,
                         type.measure = "class",
                         family = "multinomial", 
                         type.multinomial = "ungrouped")
  
  PredTrain <- predict(object = lasso.tun, 
                       newx = XTrainBalanceado, 
                       type = "class", 
                       s = "lambda.min")
  
  PredTest <- predict(object = lasso.tun, 
                      newx = XTest, 
                      type = "class", 
                      s = "lambda.min")
  
  MC.Train <- table(PredTrain, TrainBalanceado$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M1.MR_lasso <- GenerarResultadosParalelo(Metodo = "MR_lasso1", 
                                         workers = availableCores())

M2.MR_lasso <- GenerarResultadosParalelo(Metodo = "MR_lasso2", 
                                         workers = availableCores())

M3.MR_lasso <- GenerarResultadosParalelo(Metodo = "MR_lasso3", 
                                         workers = availableCores())




