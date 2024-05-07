
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

library(MASS)
library(NetDA)

library(ggplot2)

library(readr)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

source("../../Network Discriminant Analysis/NetQDA.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

QDA.1 <- function(Train, Test) {
  
  QDA <- qda(formula = Clase~.,
             data = Train)
  
  PredTrain <- predict(object = QDA,
                       newdata = Train)$class
  
  PredTest <- predict(object = QDA,
                      newdata = Test)$class
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

QDA.2 <- function(Train, Test) {
  
  X_Train <- Train %>% 
    dplyr::select(-Clase) %>% 
    as.matrix()
  
  X_Test <- Test %>% 
    dplyr::select(-Clase) %>% 
    as.matrix()
  
  Y_Train <- Train %>% 
    mutate(Clase = case_when(
      Clase == "GTEX_B" ~ 1,
      Clase == "TCGA_BLGG" ~ 2,
      Clase == "TCGA_GM" ~ 3)) %>% 
    dplyr::select(Clase) %>% 
    as.matrix()
  
  Y_Test <- Test %>% 
    mutate(Clase = case_when(
      Clase == "GTEX_B" ~ 1,
      Clase == "TCGA_BLGG" ~ 2,
      Clase == "TCGA_GM" ~ 3)) %>% 
    dplyr::select(Clase) %>% 
    as.matrix()
  
  NetQDA.Train <- NetDA(X = X_Train,
                        Y = Y_Train,
                        method = 2,
                        X_test = X_Train)
  
  NetQDA.Test <- NetDA(X = X_Train,
                       Y = Y_Train,
                       method = 2,
                       X_test = X_Test)
  
  MC.Train <- table(Y_Train, NetQDA.Train$yhat)
  MC.Test <- table(Y_Test, NetQDA.Test$yhat)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

QDA.3 <- function(Train, Test) {
  
  rho.tune <- tune_rho_NetQDA(formula = Clase~.,
                              datos = Train,
                              prior_prob = c(1/3,1/3,1/3),
                              rhos = seq(0.1, 1, by = 0.001))
  
  NetQDA <- NetQDA(formula = Clase~., 
                   datos = Train,
                   rho = rho.tune$best_rho,
                   prior_prob = c(1/3,1/3,1/3))
  
  PredTrain <- Predict.NetQDA(object = NetQDA,
                              NewData = dplyr::select(Train, -Clase))
  
  PredTest <- Predict.NetQDA(object = NetQDA,
                             NewData = dplyr::select(Test, -Clase))
  
  MC.Train <- table(Train$Clase, PredTrain$clase_predicha$.)
  MC.Test <- table(Test$Clase, PredTest$clase_predicha$.)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

# Resultados --------------------------------------------------------------

M3.QDA.1 <- Evaluacion(Metodo = "QDA.1", 
                       workers = availableCores())

M3.QDA.2 <- Evaluacion(Metodo = "QDA.2", 
                       workers = availableCores())

M3.QDA.3 <- Evaluacion(Metodo = "QDA.3", 
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.QDA.1 <- M3.QDA.1[["Global"]] %>% 
  mutate(Modelo = "M3",
         Nombre = "QDA.1")

G.QDA.2 <- M3.QDA.2[["Global"]] %>% 
  mutate(Modelo = "M3",
         Nombre = "QDA.2")

G.QDA.3 <- M3.QDA.3[["Global"]] %>% 
  mutate(Modelo = "M3",
         Nombre = "QDA.3")

M3 <- bind_rows(G.QDA.1, G.QDA.2, G.QDA.3) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M3,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.QDA.1$TestGlobal)
mean(G.QDA.2$TestGlobal)
mean(G.QDA.3$TestGlobal)

write.csv(x = M3,
          file = "Modelo3.csv",
          row.names = FALSE)
