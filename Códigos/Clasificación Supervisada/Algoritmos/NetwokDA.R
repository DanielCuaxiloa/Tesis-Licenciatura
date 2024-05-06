
#################################
# Network Discriminant Analysis #
# - Network LDA                 #
# - Network QDA                 #
#################################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(ggplot2)

library(readr)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

source("../../Network Discriminant Analysis/NetQDA.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

Network.QDA <- function(Train, Test) {
  
  NetQDA <- NetQDA(formula = Clase~., 
                   datos = Train,
                   rho = 0.02)
  
  PredTrain <- Predict.NetQDA(object = NetQDA,
                              NewData = select(Train, -Clase))
  
  PredTest <- Predict.NetQDA(object = NetQDA,
                             NewData = select(Test, -Clase))
  
  MC.Train <- table(Train$Clase, PredTrain$clase_predicha$.)
  MC.Test <- table(Test$Clase, PredTest$clase_predicha$.)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M7.Network.QDA <- Evaluacion(Metodo = "Network.QDA",
                             workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.Network.QDA <- M7.Network.QDA[["Global"]] %>% 
  mutate(Modelo = "M7",
         Nombre = "Network.QDA")

M7 <- bind_rows(G.Network.QDA) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M7,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.Network.QDA$TestGlobal)

write.csv(x = M2,
          file = "Modelo2.csv",
          row.names = FALSE)



