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

source("../../FuncionesAuxiliares.R")

source("../../../UGGM_V1/UGGM_QDA.R")


# Conjuntos Train y Test --------------------------------------------------

load("../../Folds_3Clases.RData")


# Esquemas de clasificación -----------------------------------------------


QDA.2 <- function(Train, Test) {
  
  rho.tune <- tune.rho(formula = Clase~.,
                       data = Train,
                       rhos = list(seq(0.01, 0.1, by = 0.01),
                                   seq(0.01, 0.1, by = 0.01),
                                   seq(0.01, 0.1, by = 0.01)),
                       prior = c(1/3, 1/3, 1/3),
                       nfolds = 5)
  
  NetQDA <- UGGM_QDA(formula = Clase~., 
                     data = Train,
                     rho = c(rho.tune$best.rho$Var1,
                             rho.tune$best.rho$Var2,
                             rho.tune$best.rho$Var3),
                     prior = c(1/3, 1/3, 1/3))
  
  PredTrain <- predict.UGGM_QDA(object = NetQDA,
                                newdata = dplyr::select(Train, -Clase))
  
  PredTest <- predict.UGGM_QDA(object = NetQDA,
                               newdata = dplyr::select(Test, -Clase))
  
  MC.Train <- table(Train$Clase, PredTrain$clase$yhat)
  MC.Test <- table(Test$Clase, PredTest$clase$yhat)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

M3.QDA.2 <- Evaluacion1(Metodo = "QDA.2", 
                        workers = availableCores())

# Gráficas ----------------------------------------------------------------

G.QDA.2 <- M3.QDA.2[["Global"]] %>% 
  mutate(Modelo = "Quadratic Discriminant Analysis",
         Nombre = "QDA 2")

mean(G.QDA.2$TestGlobal)

write.csv(x = M3,
          file = "Modelo3.csv",
          row.names = FALSE)

# Matriz de confusión promediada ------------------------------------------

MC.M3.QDA.1 <- M3.QDA.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M3.QDA.2 <- M3.QDA.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M3.QDA.1.PROM <- Reduce("+", MC.M3.QDA.1$MC.Test) / length(MC.M3.QDA.1$MC.Test)
MC.M3.QDA.1.PROM <- round(t(apply(MC.M3.QDA.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M3.QDA.2.PROM <- Reduce("+", MC.M3.QDA.2$MC.Test) / length(MC.M3.QDA.2$MC.Test)
MC.M3.QDA.2.PROM <- round(t(apply(MC.M3.QDA.2.PROM, 1, function(x) x / sum(x) * 100)),2)
