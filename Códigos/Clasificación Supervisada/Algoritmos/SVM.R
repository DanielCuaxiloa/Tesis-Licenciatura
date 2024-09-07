
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

SVM.1 <- function(Train, Test) {
  
  svm <- svm(formula = Clase~.,
             kernel = "radial",
             data = Train)
  
  PredTrain <- predict(object = svm, 
                       newdata = Train, 
                       type = "response")
  
  PredTest <- predict(object = svm, 
                      newdata = Test, 
                      type = "response")
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

SVM.2 <- function(Train, Test) {
  
  svm.tune <- tune(svm, Clase~., 
                   data = Train, 
                   kernel = "radial", 
                   ranges = list(cost = seq(from = 1, to = 100, by = 10), 
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
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M4.SVM.1 <- Evaluacion(Metodo = "SVM.1", 
                       workers = availableCores())

M4.SVM.2 <- Evaluacion(Metodo = "SVM.2", 
                       workers = availableCores())


# Grafica -----------------------------------------------------------------

G.SVM.1 <- M4.SVM.1[["Global"]] %>%
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 1")

G.SVM.2 <- M4.SVM.2[["Global"]] %>% 
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 2")

M4 <- bind_rows(G.SVM.1, G.SVM.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M4,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.SVM.1$TestGlobal)
mean(G.SVM.2$TestGlobal)

write.csv(x = M4,
          file = "Modelo4.csv",
          row.names = FALSE)








