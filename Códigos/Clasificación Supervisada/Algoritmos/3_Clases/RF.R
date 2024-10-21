
#################
# Random Forest #
# RF            #
#################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(ranger)

library(ggplot2)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

RF.1 <- function(Train, Test) {
  
  RF <- ranger(formula = Clase~.,
               data = Train,
               importance = "impurity",
               probability = FALSE)
  
  PredTrain <- predict(object = RF,
                       data = Train)
  
  PredTest <- predict(object = RF,
                      data = Test)
  
  MC.Train <- table(Train$Clase, PredTrain$predictions)
  MC.Test <- table(Test$Clase, PredTest$predictions)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

RF.2 <- function(Train, Test) {
  
  malla_hyper <- expand.grid(
    mtry       = seq(from = 1, to = 11, by = 1),
    node_size  = c(1,10,15),
    num.trees  = c(50,100,500)    
  )

  malla_hyper$OOBerr <- NA
  
  for(i in 1:nrow(malla_hyper)) {
    rf <- ranger(
      formula        = Clase~.,
      data           = Train,
      num.trees      = malla_hyper$num.trees[i],
      mtry           = malla_hyper$mtry[i],
      min.node.size  = malla_hyper$node_size[i],
      importance = 'impurity')
    malla_hyper$OOBerr[i] <- rf$prediction.error
  }
  
  position <- which.min(malla_hyper$OOBerr) 
  
  RF.Tune <- ranger(Clase~.,
                    data = Train, 
                    num.trees = malla_hyper$num.trees[position],
                    min.node.size = malla_hyper$node_size[position], 
                    mtry = malla_hyper$mtry[position],
                    importance = 'impurity', 
                    probability = FALSE)
  
  PredTrain <- predict(object = RF.Tune,
                       data = Train)
  
  PredTest <- predict(object = RF.Tune,
                      data = Test)
  
  MC.Train <- table(PredTrain$predictions, Train$Clase)
  MC.Test <- table(PredTest$predictions, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
}


# Resultados --------------------------------------------------------------

M5.RF.1 <- Evaluacion(Metodo = "RF.1",
                      workers = availableCores())

M5.RF.2 <- Evaluacion(Metodo = "RF.2",
                      workers = availableCores())


# Grafica -----------------------------------------------------------------

G.RF.1 <- M5.RF.1[["Global"]] %>%
  mutate(Modelo = "Random Forest",
         Nombre = "RF 1")

G.RF.2 <- M5.RF.2[["Global"]] %>%
  mutate(Modelo = "Random Forest",
         Nombre = "RF 2")

M5 <- bind_rows(G.RF.1, G.RF.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M5,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.RF.1$TestGlobal)
mean(G.RF.2$TestGlobal)

write.csv(x = M5,
          file = "Modelo5.csv",
          row.names = FALSE)




