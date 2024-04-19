
##########################
# Support Vector Machine #
# SVM                    #
##########################


library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(e1071)

library(ggplot2)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

M4.SVM_1 <- function(Train, Test) {
  
  svm <- svm(formula = Clase~.,
             kernel = "radial",
             data = Train)
  
  PredTrain <- predict(object = svm, 
                       newdata = Train, 
                       type = "response")
  
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

M4.SVM_2 <- function(Train, Test) {
  
  svm.tune <- tune(svm, Clase~., 
                   data = Train, 
                   kernel = "radial", 
                   ranges = list(cost = seq(from = 1, to = 5, by = 1), 
                                 gamma = seq(from = 0.1, to = 1, by = 0.1)))
  
  svm <- svm(formula = Clase~.,
             cost = svm.tune$best.parameters[[1]], 
             gamma = svm.tune$best.parameters[[2]],
             data = Train)
  
  PredTrain <- predict(object = svm, 
                       newdata = Train, 
                       type = "response")
  
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
  
  MC.Train <- table(PredTrain, Train$Clase)
  MC.Test <- table(PredTest, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


M4.SVM_1 <- Evaluacion(Metodo = "M4.SVM_1", 
                       workers = availableCores())

M4.SVM_2 <- Evaluacion(Metodo = "M4.SVM_2", 
                       workers = availableCores())



# Grafica -----------------------------------------------------------------

G.SVM_1 <- M4.SVM_1[["Global"]] %>%
  mutate(Modelo = "M4",
         Nombre = "SVM_1")

G.SVM_2 <- M4.SVM_2[["Global"]] %>% 
  mutate(Modelo = "M4",
         Nombre = "SVM_2")

M4 <- bind_rows(G.SVM_1, G.SVM_2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M4,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.SVM_1$TestGlobal)
mean(G.SVM_2$TestGlobal)

write.csv(x = M4,
          file = "Modelo_4.csv",
          row.names = FALSE)








