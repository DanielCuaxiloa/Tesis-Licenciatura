
###################################
# Multinomial Logistic Regression #
# MLR Elastic-Net                 #
###################################


library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)

library(rsample)

library(purrr)
library(furrr)

library(VGAM)
library(glmnet)

library(ggplot2)

library(readr)

library(tictoc)


# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Algoritmos de clasificación ---------------------------------------------

MLR.1 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~., 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~., 
                        data = Test)[,-1]
  
  MLR <- glmnet(x = XTrain, 
                y = YTrain,
                family = "multinomial", 
                type.multinomial = "ungrouped",
                lambda = 0)
  
  PredTrain <- predict(object = MLR, 
                       newx = XTrain, 
                       type = "class")
  
  PredTest <- predict(object = MLR, 
                      newx = XTest, 
                      type = "class")
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

MLR.2 <- function(Train, Test) {  
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  alphas <- seq(0, 1, by = 0.5)
  
  mejor_modelo <- NULL
  mejor_alpha <- NULL
  mejor_lambda <- NULL
  mejor_error <- Inf
  
  for (a in alphas) {
    modelo_tuneo <- cv.glmnet(
      x = XTrain, 
      y = YTrain, 
      nfolds = 3,
      alpha = a,
      type.measure = "class",
      family = "multinomial", 
      type.multinomial = "ungrouped"
    )
    
    error_actual <- min(modelo_tuneo$cvm)
    lambda_actual <- modelo_tuneo$lambda.min
    
    if (error_actual < mejor_error) {
      mejor_error <- error_actual
      mejor_modelo <- modelo_tuneo
      mejor_alpha <- a
      mejor_lambda <- lambda_actual
    }
  }
  
  PredTrain <- predict(object = mejor_modelo, 
                       newx = XTrain, 
                       type = "class", 
                       s = mejor_lambda)
  
  PredTest <- predict(object = mejor_modelo, 
                      newx = XTest, 
                      type = "class", 
                      s = mejor_lambda)
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
}


# Resultados --------------------------------------------------------------

tic()
M1.MLR.1 <- Evaluacion1(Metodo = "MLR.1", 
                        workers = availableCores())
TE_MLR.1 <- toc(log = TRUE)

tic()
M1.MLR.2 <- Evaluacion1(Metodo = "MLR.2",
                        workers = availableCores())
TE_MLR.2 <- toc(log = TRUE)


# Gráficas ----------------------------------------------------------------

G.MLR.1 <- M1.MLR.1[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 1")

G.MLR.2 <- M1.MLR.2[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 2")

M1 <- bind_rows(G.MLR.1, G.MLR.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

G.TE_MLR.1 <- data.frame(TE = TE_MLR.1[["callback_msg"]]) %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 1")

G.TE_MLR.2 <- data.frame(TE = TE_MLR.2[["callback_msg"]]) %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 2")

TE_M2 <- bind_rows(G.TE_MLR.1, G.TE_MLR.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

write.csv(x = M1,
          file = "Modelo1.csv",
          row.names = FALSE)

write.csv(x = TE_M2,
          file = "TE_Modelo1.csv",
          row.names = FALSE)


# Matriz de confusión promediada ------------------------------------------

MC.M1.MLR.1 <- M1.MLR.1[["MatricesConfusion"]] %>% 
  transpose()

MC.M1.MLR.2 <- M1.MLR.2[["MatricesConfusion"]] %>% 
  transpose()

MC.M1.MLR.1.PROM <- Reduce("+", MC.M1.MLR.1$MC.Test) / length(MC.M1.MLR.1$MC.Test)
MC.M1.MLR.1.PROM <- round(t(apply(MC.M1.MLR.1.PROM, 1, function(x) x / sum(x) * 100)),2)

MC.M1.MLR.2.PROM <- Reduce("+", MC.M1.MLR.2$MC.Test) / length(MC.M1.MLR.2$MC.Test)
MC.M1.MLR.2.PROM <- round(t(apply(MC.M1.MLR.2.PROM, 1, function(x) x / sum(x) * 100)),2)


# MC Mapa de Calor --------------------------------------------------------

MC.M1.MLR.1.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(88.41, 10.58, 1.01),
  TCGA_LGG = c(4.49, 89.25, 6.25),
  TCGA_GM = c(4.04, 29.97, 65.99)) %>%
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
         Modelo = as.factor("Multinomial Logistic Regression"),
         Nombre = as.factor("Modelo 1")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

MC.M1.MLR.2.DF <- tibble(
  PredTest = c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"),
  GTEX_Brain = c(89.13, 10.09, 0.78),
  TCGA_LGG = c(4.15, 89.30, 6.55),
  TCGA_GM = c(3.34, 29.34, 67.32)) %>%
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
         Modelo = as.factor("Multinomial Logistic Regression"),
         Nombre = as.factor("Modelo 2")) %>% 
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

M1.MC <- bind_rows(MC.M1.MLR.1.DF, MC.M1.MLR.2.DF)

write.csv(x = M1.MC,
          file = "MC1.csv",
          row.names = FALSE)
