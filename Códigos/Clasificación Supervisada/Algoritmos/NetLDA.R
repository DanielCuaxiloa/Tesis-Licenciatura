
########################################
# Network Linear Discriminant Analysis #
# NetLDA                               #
########################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(NetDA)

library(ggplot2)

library(readr)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

NetLDA <- function(Train, Test) {
  
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

NetQDA <- function(Train, Test) {
  
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


# Resultados --------------------------------------------------------------

M6.NetLDA <- Evaluacion(Metodo = "NetLDA",
                        workers = availableCores())

M6.NetQDA <- Evaluacion(Metodo = "NetQDA",
                        workers = availableCores())

# Gráficas ----------------------------------------------------------------

G.NetLDA <- M6.NetLDA[["Global"]] %>% 
  mutate(Modelo = "M6",
         Nombre = "NetLDA")

G.NetQDA <- M6.NetQDA[["Global"]] %>% 
  mutate(Modelo = "M6",
         Nombre = "NetQDA")

M6 <- bind_rows(G.NetLDA, G.NetQDA) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M6,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.NetLDA$TestGlobal)
mean(G.NetQDA$TestGlobal)

write.csv(x = M6,
          file = "Modelo_6.csv",
          row.names = FALSE)



