
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

library(MASS)
library(NetDA)

library(ggplot2)

library(readr)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

LDA.1 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~., 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class

  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LDA.2 <- function(Train, Test) {
  
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
  
  NetLDA.Train <- NetDA(X = X_Train,
                        Y = Y_Train,
                        method = 1,
                        X_test = X_Train)
  
  NetLDA.Test <- NetDA(X = X_Train,
                       Y = Y_Train,
                       method = 1,
                       X_test = X_Test)
  
  MC.Train <- table(Y_Train, NetLDA.Train$yhat)
  MC.Test <- table(Y_Test, NetLDA.Test$yhat)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M2.LDA.1 <- Evaluacion(Metodo = "LDA.1",
                       workers = availableCores())

M2.LDA.2 <- Evaluacion(Metodo = "LDA.2", 
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.LDA.1 <- M2.LDA.1[["Global"]] %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 1")

G.LDA.2 <- M2.LDA.2[["Global"]] %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 2")

M2 <- bind_rows(G.LDA.1, G.LDA.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M2,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.LDA.1$TestGlobal)
mean(G.LDA.2$TestGlobal)

write.csv(x = M2,
          file = "Modelo2.csv",
          row.names = FALSE)


# Matriz de confusión promediada ------------------------------------------

MC.M2.LDA.1 <- M2.LDA.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M2.LDA.2 <- M2.LDA.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M2.LDA.1.PROM <- Reduce("+", MC.M2.LDA.1$MC.Test) / length(MC.M2.LDA.1$MC.Test)
MC.M2.LDA.1.PROM <- round(t(apply(MC.M2.LDA.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M2.LDA.2.PROM <- Reduce("+", MC.M2.LDA.2$MC.Test) / length(MC.M2.LDA.2$MC.Test)
MC.M2.LDA.2.PROM <- round(t(apply(MC.M2.LDA.2.PROM, 1, function(x) x / sum(x) * 100)),2)






