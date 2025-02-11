
###################################
# 3 Clases                        #
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


# Funciones Auxiliares ----------------------------------------------------

source("../../FuncionesAuxiliares.R")


# Conjuntos Train y Test --------------------------------------------------

load("../../Folds_3Clases.RData")


# Algoritmos de clasificacion ---------------------------------------------

MLR.1 <- function(Train, Test) {
  
  lrm <- vglm(formula = Clase~.,
              family = multinomial(),
              data = Train)
  
  PredTrain_Probs <- predict(object = lrm,
                             newdata = select(Train,-Clase),
                             type = "response")
  
  PredTrain_Class <- data.frame(Clase = levels(Train$Clase)[max.col(PredTrain_Probs)])
  
  PredTest_Probs <- predict(object = lrm,
                            newdata = select(Test,-Clase),
                            type = "response")
  
  PredTest_Class <- data.frame(Clase = levels(Test$Clase)[max.col(PredTest_Probs)])
  
  MC.Train <- table(Train$Clase, PredTrain_Class$Clase)
  MC.Test <- table(Test$Clase, PredTest_Class$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

MLR.2 <- function(Train, Test) {
  
  XTrain <- model.matrix(Clase~.^2, 
                         data = Train)[,-1]
  YTrain <- Train$Clase
  
  XTest <- model.matrix(Clase~.^2, 
                        data = Test)[,-1]
  
  lasso.tun <- cv.glmnet(x = XTrain, 
                         y = YTrain, 
                         nfolds = 5,
                         alpha = 1,
                         #lambda = seq(from = 0, to = 20, by = 0.5),
                         type.measure = "class",
                         family = "multinomial", 
                         type.multinomial = "ungrouped")
  
  PredTrain <- predict(object = lasso.tun, 
                       newx = XTrain, 
                       type = "class", 
                       s = "lambda.min")
  
  PredTest <- predict(object = lasso.tun, 
                      newx = XTest, 
                      type = "class", 
                      s = "lambda.min")
  
  MC.Train <- table(Train$Clase, PredTrain)
  MC.Test <- table(Test$Clase, PredTest)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M1.MLR.1 <- Evaluacion1(Metodo = "MLR.1", 
                        workers = availableCores())

M1.MLR.2 <- Evaluacion1(Metodo = "MLR.2",
                        workers = availableCores())


# Graficas ----------------------------------------------------------------

G.MLR.1 <- M1.MLR.1[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 1")

G.MLR.2 <- M1.MLR.2[["Global"]] %>% 
  mutate(Modelo = "Multinomial Logistic Regression",
         Nombre = "MLR 2")

M1 <- bind_rows(G.MLR.1, G.MLR.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

write.csv(x = M1,
          file = "Modelo1.csv",
          row.names = FALSE)


# Matriz de confusion promediada ------------------------------------------

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
