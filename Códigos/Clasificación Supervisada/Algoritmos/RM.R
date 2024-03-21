
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

library(VGAM)
library(glmnet)

library(ggplot2)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Algoritmos de clasificación ---------------------------------------------

LRM_1 <- function(Train, Test) {
  
  lrm <- vglm(formula = Clase~.,
              family = multinomial(),
              data = Train)
  
  PredTrain_Probs <- predict(object = lrm,
                             newdata = select(Train,-Clase),
                             type = "response")
  
  PredTrain_Class <- data.frame(Clase = levels(Train$Clase)[max.col(PredTrain_Probs)])
  
  PredTest_Probs <- predict(object = lrm,
                            newdata = select(Test,-Clase),
                            type = "response")
  
  PredTest_Class <- data.frame(Clase = levels(Test$Clase)[max.col(PredTest_Probs)])
  
  MC.Train <- table(PredTrain_Class$Clase, Train$Clase)
  MC.Test <- table(PredTest_Class$Clase, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LRM_2 <- function(Train, Test) {
  
  lrm <- vglm(formula = Clase~.^2,
              family = multinomial(),
              data = Train)
  
  PredTrain_Probs <- predict(object = lrm,
                             newdata = select(Train,-Clase),
                             type = "response")
  
  PredTrain_Class <- data.frame(Clase = levels(Train$Clase)[max.col(PredTrain_Probs)])
  
  PredTest_Probs <- predict(object = lrm,
                            newdata = select(Test,-Clase),
                            type = "response")
  
  PredTest_Class <- data.frame(Clase = levels(Test$Clase)[max.col(PredTest_Probs)])
  
  MC.Train <- table(PredTrain_Class$Clase, Train$Clase)
  MC.Test <- table(PredTest_Class$Clase, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LRM_lasso_1 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~., 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~., 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 10,
                         alpha = 0.2,
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

LRM_lasso_2 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 10,
                         alpha = 0.2,
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

LRM_lasso_SMOTE <- function(Train, Test) {
  
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

M1.LRM_1 <- GenerarResultadosParalelo(Metodo = "LRM_1", 
                                    workers = availableCores())

M1.LRM_2 <- GenerarResultadosParalelo(Metodo = "LRM_2", 
                                      workers = availableCores())

M1.LRM_lasso_1 <- GenerarResultadosParalelo(Metodo = "LRM_lasso_1", 
                                            workers = availableCores())

M1.LRM_lasso_2 <- GenerarResultadosParalelo(Metodo = "LRM_lasso_2", 
                                            workers = availableCores())


# Gráficas ----------------------------------------------------------------

LRM_1 <- M1.LRM_1[["Global"]] %>% 
  mutate(Modelo = "M1.LRM_1")

LRM_2 <- M1.LRM_2[["Global"]] %>% 
  mutate(Modelo = "M1.LRM_2")

LRM_lasso_1 <- M1.LRM_lasso_1[["Global"]] %>% 
  mutate(Modelo = "M1.LRM_lasso_1")

LRM_lasso_2 <- M1.LRM_lasso_2[["Global"]] %>% 
  mutate(Modelo = "M1.LRM_lasso_2")

M1 <- bind_rows(LRM_1, LRM_2, LRM_lasso_1, LRM_lasso_2) %>% 
  mutate(Modelo = as.factor(Modelo))

ggplot(data = M1,
       mapping = aes(x = Modelo, y = TestGlobal_KCV)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

