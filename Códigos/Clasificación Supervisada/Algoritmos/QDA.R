
###################################
# Quadratic Discriminant Analysis #
# QDA                             #
###################################


library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)

library(rsample)

library(purrr)
library(furrr)

library(MASS)

library(ggplot2)

library(readr)

library(tictoc)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

source("../../UGGM/Modelo_UGGM-QDA.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

QDA.1 <- function(Train, Test) {
  
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

QDA.2 <- function(Train, Test) {
  
  rho.tune <- tune.rho(formula = Clase~.,
                       data = Train,
                       rhos = seq(0.1, 1, by = 0.01),
                       #prior = c(1/3, 1/3, 1/3),
                       nfolds = 5)
  
  NetQDA <- UGGM_QDA(formula = Clase~., 
                     data = Train,
                     rho = rho.tune$best.rho,
                     prior = c(1/3, 1/3, 1/3))
  
  PredTrain <- predict.UGGM_QDA(object = NetQDA,
                                newdata = dplyr::select(Train, -Clase))
  
  PredTest <- predict.UGGM_QDA(object = NetQDA,
                               newdata = dplyr::select(Test, -Clase))
  
  MC.Train <- table(Train$Clase, PredTrain$clase$yhat)
  MC.Test <- table(Test$Clase, PredTest$clase$yhat)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

tic()
M3.QDA.1 <- Evaluacion1(Metodo = "QDA.1", 
                        workers = availableCores())
TE_QDA.1 <- toc(log = TRUE)

tic()
M3.QDA.2 <- Evaluacion1(Metodo = "QDA.2", 
                        workers = availableCores())
TE_QDA.2 <- toc(log = TRUE)


# Gráficas ----------------------------------------------------------------

G.QDA.1 <- M3.QDA.1[["Global"]] %>% 
  mutate(Modelo = "Quadratic Discriminant Analysis",
         Nombre = "QDA 1")

G.QDA.2 <- M3.QDA.2[["Global"]] %>% 
  mutate(Modelo = "Quadratic Discriminant Analysis",
         Nombre = "QDA 2")

M3 <- bind_rows(G.QDA.1, G.QDA.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

G.TE_QDA.1 <- data.frame(TE = TE_QDA.1[["callback_msg"]]) %>% 
  mutate(Modelo = "Quadratic Discriminant Analysis",
         Nombre = "QDA 1")

G.TE_QDA.2 <- data.frame(TE = TE_QDA.2[["callback_msg"]]) %>% 
  mutate(Modelo = "Quadratic Discriminant Analysis",
         Nombre = "QDA 2")

TE_M3 <- bind_rows(G.TE_QDA.1, G.TE_QDA.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

write.csv(x = M3,
          file = "Modelo3.csv",
          row.names = FALSE)

write.csv(x = TE_M3,
          file = "TE_Modelo3.csv",
          row.names = FALSE)


# Matriz de confusión promediada ------------------------------------------

MC.M3.QDA.1 <- M3.QDA.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M3.QDA.2 <- M3.QDA.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M3.QDA.1.PROM <- Reduce("+", MC.M3.QDA.1$MC.Test) / length(MC.M3.QDA.1$MC.Test)
MC.M3.QDA.1.PROM <- round(t(apply(MC.M3.QDA.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M3.QDA.2.PROM <- Reduce("+", MC.M3.QDA.2$MC.Test) / length(MC.M3.QDA.2$MC.Test)
MC.M3.QDA.2.PROM <- round(t(apply(MC.M3.QDA.2.PROM, 1, function(x) x / sum(x) * 100)),2)


# MC Mapa de Calor --------------------------------------------------------

MC.M3.QDA.1.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(89.88, 8.57, 1.55),
  TCGA_LGG = c(5.41, 87.34, 7.25),
  TCGA_GM = c(2.14, 25.27, 72.59)) %>%
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
         Modelo = as.factor("Quadratic Discriminant Analysis"),
         Nombre = as.factor("Modelo 1")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

MC.M3.QDA.2.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(91.06, 7.23, 1.71),
  TCGA_LGG = c(9.68, 77.55, 12.76),
  TCGA_GM = c(3.49, 13.37, 83.13)) %>%
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
         Modelo = as.factor("Quadratic Discriminant Analysis"),
         Nombre = as.factor("Modelo 2")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

M3.MC <- bind_rows(MC.M3.QDA.1.DF, MC.M3.QDA.2.DF)

write.csv(x = M3.MC,
          file = "MC3.csv",
          row.names = FALSE)
