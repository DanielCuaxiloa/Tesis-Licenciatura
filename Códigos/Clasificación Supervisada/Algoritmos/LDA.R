
################################
# Linear Discriminant Analysis #
# LDA                          #
################################


library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)

library(rsample)

library(purrr)
library(furrr)

library(MASS)
library(NetDA)

library(ggplot2)

library(readr)

library(tictoc)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

LDA.1 <- function(Train, Test) {
  
  LDA <- lda(formula = Clase~., 
             data = Train)
  
  PredTrain <- predict(object = LDA,
                       newdata = Train)$class

  PredTest <- predict(object = LDA,
                      newdata = Test)$class
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

LDA.2 <- function(Train, Test) {
  
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


# Resultados --------------------------------------------------------------

tic()
M2.LDA.1 <- Evaluacion1(Metodo = "LDA.1",
                        workers = availableCores())
TE_LDA.1 <- toc(log = TRUE)

tic()
M2.LDA.2 <- Evaluacion1(Metodo = "LDA.2", 
                        workers = availableCores())
TE_LDA.2 <- toc(log = TRUE)


# Gráficas ----------------------------------------------------------------

G.LDA.1 <- M2.LDA.1[["Global"]] %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 1")

G.LDA.2 <- M2.LDA.2[["Global"]] %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 2")

M2 <- bind_rows(G.LDA.1, G.LDA.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

G.TE_LDA.1 <- data.frame(TE = TE_LDA.1[["callback_msg"]]) %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 1")

G.TE_LDA.2 <- data.frame(TE = TE_LDA.2[["callback_msg"]]) %>% 
  mutate(Modelo = "Linear Discriminant Analysis",
         Nombre = "LDA 2")

TE_M2 <- bind_rows(G.TE_LDA.1, G.TE_LDA.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

write.csv(x = M2,
          file = "Modelo2.csv",
          row.names = FALSE)

write.csv(x = TE_M2,
          file = "TE_Modelo2.csv",
          row.names = FALSE)


# Promedios matriz de confusión -------------------------------------------

MC.M2.LDA.1 <- M2.LDA.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M2.LDA.2 <- M2.LDA.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M2.LDA.1.PROM <- Reduce("+", MC.M2.LDA.1$MC.Test) / length(MC.M2.LDA.1$MC.Test)
MC.M2.LDA.1.PROM <- round(t(apply(MC.M2.LDA.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M2.LDA.2.PROM <- Reduce("+", MC.M2.LDA.2$MC.Test) / length(MC.M2.LDA.2$MC.Test)
MC.M2.LDA.2.PROM <- round(t(apply(MC.M2.LDA.2.PROM, 1, function(x) x / sum(x) * 100)),2)


# MC Mapa de Calor --------------------------------------------------------

MC.M2.LDA.1.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(87.56, 11.70, 0.74),
  TCGA_LGG = c(6.37, 85.52, 8.12),
  TCGA_GM = c(2.05, 27.20, 70.75)) %>%
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
         Modelo = as.factor("Linear Discriminant Analysis"),
         Nombre = as.factor("Modelo 1")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

MC.M2.LDA.2.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(81.77, 17.53, 0.71),
  TCGA_LGG = c(3.93, 91.20, 4.88),
  TCGA_GM = c(2.20, 40.66, 57.14)) %>%
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
         Modelo = as.factor("Linear Discriminant Analysis"),
         Nombre = as.factor("Modelo 2")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

M2.MC <- bind_rows(MC.M2.LDA.1.DF, MC.M2.LDA.2.DF)

write.csv(x = M2.MC,
          file = "MC2.csv",
          row.names = FALSE)
