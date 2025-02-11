
#################
# 3 Clases      #
# Random Forest #
# RF            #
#################


library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)

library(rsample)

library(purrr)
library(furrr)

library(ranger)

library(ggplot2)

# Funciones Auxiliares ----------------------------------------------------

source("../../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../../Folds_3Clases.RData")


# Esquemas de clasificacion -----------------------------------------------

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

M5.RF.1 <- Evaluacion1(Metodo = "RF.1",
                       workers = availableCores())

M5.RF.2 <- Evaluacion1(Metodo = "RF.2",
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

write.csv(x = M5,
          file = "Modelo5.csv",
          row.names = FALSE)


# Matriz de confusion promediada ------------------------------------------

MC.M5.RF.1 <- M5.RF.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M5.RF.2 <- M5.RF.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M5.RF.1.PROM <- Reduce("+", MC.M5.RF.1$MC.Test) / length(MC.M5.RF.1$MC.Test)
MC.M5.RF.1.PROM <- round(t(apply(MC.M5.RF.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M5.RF.2.PROM <- Reduce("+", MC.M5.RF.2$MC.Test) / length(MC.M5.RF.2$MC.Test)
MC.M5.RF.2.PROM <- round(t(apply(MC.M5.RF.2.PROM, 1, function(x) x / sum(x) * 100)),2)


# MC Mapa de Calor --------------------------------------------------------

MC.M5.RF.1.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(87.26, 11.78, 0.95),
  TCGA_LGG = c(4.30, 89.89, 5.80),
  TCGA_GM = c(2.02, 31.20, 66.78)) %>%
  pivot_longer(-PredTest, names_to = "Referencia", values_to = "Porcentaje") %>% 
  rename(Prediccion = PredTest) %>%
  mutate(Referencia = as.factor(Referencia),
         Prediccion = as.factor(Prediccion)) %>% 
  mutate(Prediccion = recode(Prediccion, 
                             "GTEX_Brain" = "GTEX Brain", 
                             "TCGA_LGG" = "TCGA LGG", 
                             "TCGA_GM" = "TCGA GM"),
         Referencia = recode(Referencia,
                             "GTEX_Brain" = "GTEX Brain",
                             "TCGA_LGG" = "TCGA LGG",
                             "TCGA_GM" = "TCGA GM"),
         Modelo = as.factor("Random Forest"),
         Nombre = as.factor("Modelo 1")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

MC.M5.RF.2.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(90.77, 8.14, 1.09),
  TCGA_LGG = c(6.22, 84.31, 9.47),
  TCGA_GM = c(1.96, 20.96, 77.08)) %>%
  pivot_longer(-PredTest, names_to = "Referencia", values_to = "Porcentaje") %>% 
  rename(Prediccion = PredTest) %>%
  mutate(Referencia = as.factor(Referencia),
         Prediccion = as.factor(Prediccion)) %>% 
  mutate(Prediccion = recode(Prediccion, 
                             "GTEX_Brain" = "GTEX Brain", 
                             "TCGA_LGG" = "TCGA LGG", 
                             "TCGA_GM" = "TCGA GM"),
         Referencia = recode(Referencia,
                             "GTEX_Brain" = "GTEX Brain",
                             "TCGA_LGG" = "TCGA LGG",
                             "TCGA_GM" = "TCGA GM"),
         Modelo = as.factor("Random Forest"),
         Nombre = as.factor("Modelo 2")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

M5.MC <- bind_rows(MC.M5.RF.1.DF, MC.M5.RF.2.DF)

write.csv(x = M5.MC,
          file = "MC5.csv",
          row.names = FALSE)

