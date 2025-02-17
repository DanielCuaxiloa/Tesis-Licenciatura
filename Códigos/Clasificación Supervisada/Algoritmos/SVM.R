
##########################
# Support Vector Machine #
# SVM                    #
##########################


library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)

library(rsample)

library(purrr)
library(furrr)

library(e1071)

library(ggplot2)

library(tictoc)


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
                   tunecontrol = tune.control(cross = 5),
                   ranges = list(cost = 2^(2:6), gamma = 2^(-6:-2)))
  
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

tic()
M4.SVM.1 <- Evaluacion1(Metodo = "SVM.1", 
                        workers = availableCores())
TE_SVM.1 <- toc(log = TRUE)

tic()
M4.SVM.2 <- Evaluacion1(Metodo = "SVM.2", 
                        workers = availableCores())
TE_SVM.2 <- toc(log = TRUE)


# Gráfica -----------------------------------------------------------------

G.SVM.1 <- M4.SVM.1[["Global"]] %>%
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 1")

G.SVM.2 <- M4.SVM.2[["Global"]] %>% 
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 2")

M4 <- bind_rows(G.SVM.1, G.SVM.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

G.TE_SVM.1 <- data.frame(TE = TE_SVM.1[["callback_msg"]]) %>% 
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 1")

G.TE_SVM.2 <- data.frame(TE = TE_SVM.2[["callback_msg"]]) %>% 
  mutate(Modelo = "Support Vector Machine",
         Nombre = "SVM 2")

TE_M4 <- bind_rows(G.TE_SVM.1, G.TE_SVM.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

write.csv(x = M4,
          file = "Modelo4.csv",
          row.names = FALSE)

write.csv(x = TE_M4,
          file = "TE_Modelo4.csv",
          row.names = FALSE)


# Matriz de confusión promediada ------------------------------------------

MC.M4.SVM.1 <- M4.SVM.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M4.SVM.2 <- M4.SVM.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M4.SVM.1.PROM <- Reduce("+", MC.M4.SVM.1$MC.Test) / length(MC.M4.SVM.1$MC.Test)
MC.M4.SVM.1.PROM <- round(t(apply(MC.M4.SVM.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M4.SVM.2.PROM <- Reduce("+", MC.M4.SVM.2$MC.Test) / length(MC.M4.SVM.2$MC.Test)
MC.M4.SVM.2.PROM <- round(t(apply(MC.M4.SVM.2.PROM, 1, function(x) x / sum(x) * 100)),2)


# MC Mapa de Calor --------------------------------------------------------

MC.M4.SVM.1.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(89.17, 9.95, 0.88),
  TCGA_LGG = c(3.54, 90.90, 5.56),
  TCGA_GM = c(2.41, 31.48, 66.11)) %>%
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
         Modelo = as.factor("Support Vector Machine"),
         Nombre = as.factor("Modelo 1")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

MC.M4.SVM.2.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(88.76, 10.32, 0.92),
  TCGA_LGG = c(3.39, 90.96, 5.65),
  TCGA_GM = c(2.47, 32.29, 65.24)) %>%
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
         Modelo = as.factor("Support Vector Machine"),
         Nombre = as.factor("Modelo 2")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

M4.MC <- bind_rows(MC.M4.SVM.1.DF, MC.M4.SVM.2.DF)

write.csv(x = M4.MC,
          file = "MC4.csv",
          row.names = FALSE)
