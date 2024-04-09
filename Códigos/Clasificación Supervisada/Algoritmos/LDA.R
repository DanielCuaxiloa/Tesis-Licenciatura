
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

library(ggplot2)

library(readr)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

LDA_1 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~., 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class

  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LDA_2 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~.^2, 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class
  
  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M2.LDA_1 <- Evaluacion(Metodo = "LDA_1",
                       workers = availableCores())

M2.LDA_2 <- Evaluacion(Metodo = "LDA_2", 
                       workers = availableCores())


# Gráficas ----------------------------------------------------------------

G.LDA_1 <- M2.LDA_1[["Global"]] %>% 
  mutate(Modelo = "M2",
         Nombre = "LDA_1")

G.LDA_2 <- M2.LDA_2[["Global"]] %>% 
  mutate(Modelo = "M2",
         Nombre = "LDA_2")

M2 <- bind_rows(G.LDA_1, G.LDA_2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M2,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.LDA_1$TestGlobal)
mean(G.LDA_2$TestGlobal)

write.csv(x = M2,
          file = "Modelo_2.csv",
          row.names = FALSE)

