
###################################
# Multinomial Logistic Regression #
# MLR Elastic-Net                 #
###################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(VGAM)
library(glmnet)

library(ggplot2)

library(readr)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Algoritmos de clasificación ---------------------------------------------

MLR_1 <- function(Train, Test) {
  
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

MLR_2 <- function(Train, Test) {
  
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

MLR_3 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~., 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~., 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 5,
                         alpha = 0.8,
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

MLR_4 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 5,
                         alpha = 0.8,
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


# Resultados --------------------------------------------------------------

M1.MLR_1 <- Evaluacion(Metodo = "MLR_1", 
                       workers = availableCores())

M1.MLR_2 <- Evaluacion(Metodo = "MLR_2",
                       workers = availableCores())

M1.MLR_3 <- Evaluacion(Metodo = "MLR_3",
                       workers = availableCores())

M1.MLR_4 <- Evaluacion(Metodo = "MLR_4",
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.MLR_1 <- M1.MLR_1[["Global"]] %>% 
  mutate(Modelo = "M1",
         Nombre = "MLR_1")

G.MLR_2 <- M1.MLR_2[["Global"]] %>% 
  mutate(Modelo = "M1",
         Nombre = "MLR_2")

G.MLR_3 <- M1.MLR_3[["Global"]] %>% 
  mutate(Modelo = "M1",
         Nombre = "MLR_3")

G.MLR_4 <- M1.MLR_4[["Global"]] %>% 
  mutate(Modelo = "M1",
         Nombre = "MLR_4")

M1 <- bind_rows(G.MLR_1, G.MLR_2, G.MLR_3, G.MLR_4) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M1,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.MLR_1$TestGlobal)
mean(G.MLR_3$TestGlobal)

write.csv(x = M1,
          file = "Modelo_1.csv",
          row.names = FALSE)

