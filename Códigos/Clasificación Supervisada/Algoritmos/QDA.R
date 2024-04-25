
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

library(ggplot2)

library(readr)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

QDA_1 <- function(Train, Test) {
  
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

QDA_2 <- function(Train, Test) {
  
  QDA <- qda(formula = Clase~.^2, 
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


# Resultados --------------------------------------------------------------

M3.QDA_1 <- Evaluacion(Metodo = "QDA_1", 
                       workers = availableCores())

M3.QDA_2 <- Evaluacion(Metodo = "QDA_2", 
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.QDA_1 <- M3.QDA_1[["Global"]] %>% 
  mutate(Modelo = "M3",
         Nombre = "QDA_1")

G.QDA_2 <- M3.QDA_2[["Global"]] %>% 
  mutate(Modelo = "M3",
         Nombre = "QDA_2")

M3 <- bind_rows(G.QDA_1, G.QDA_2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M3,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.QDA_1$TestGlobal)
mean(G.QDA_2$TestGlobal)

write.csv(x = M3,
          file = "Modelo_3.csv",
          row.names = FALSE)
