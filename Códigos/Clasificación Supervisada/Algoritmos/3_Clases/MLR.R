
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

MLR.1 <- function(Train, Test) {
  
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
  
  MC.Train <- table(Train$Clase, PredTrain_Class$Clase)
  MC.Test <- table(Test$Clase, PredTest_Class$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

MLR.2 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 5,
                         alpha = 1,
                         #lambda = seq(from = 0, to = 20, by = 0.5),
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
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M1.MLR.1 <- Evaluacion(Metodo = "MLR.1", 
                       workers = availableCores())

M1.MLR.2 <- Evaluacion(Metodo = "MLR.2",
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.MLR.1 <- M1.MLR.1[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 1")

G.MLR.2 <- M1.MLR.2[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 2")

M1 <- bind_rows(G.MLR.1, G.MLR.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M1,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.MLR.1$TestGlobal)
mean(G.MLR.2$TestGlobal)

write.csv(x = M1,
          file = "Modelo1.csv",
          row.names = FALSE)

