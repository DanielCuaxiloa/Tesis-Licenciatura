####################
# Neuronal Network #
# NN               #
####################

library(tibble)
library(plyr)
library(dplyr)

library(rsample)

library(purrr)
library(furrr)

library(smotefamily)

library(keras)
library(tensorflow)
library(kerastuneR)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

NN1 <- function(Train, Test) {
  
  features_Train <- Train %>% 
    select(-Clase) %>% 
    data.frame()
  
  features_Test <- Test %>% 
    select(-Clase) %>% 
    data.frame()
  
  labels_Train <- Train %>% 
    select(Clase) %>% 
    data.frame()
  
  labels_onehot <- to_categorical(y = as.numeric(unlist(labels_Train))-1)
  
  modelo <- keras_model_sequential()
  
  modelo %>%
    layer_dense(units = 20, 
                activation = "relu", 
                input_shape = dim(features_Train)[2]) %>%
    layer_dense(units = 20,
                activation = "relu") %>% 
    layer_dense(units = 3, 
                activation = "softmax") %>% 
    compile(loss = "categorical_crossentropy",
                     optimizer = "adam",
                     metrics = c("accuracy"))
  
  historia <- modelo %>% 
    fit(x = as.matrix(features_Train), 
        y = labels_onehot,
        epochs = 100,
        batch_size = 32)
  
  y_Train <- modelo %>% 
    predict(as.matrix(features_Train)) %>% 
    k_argmax() %>% 
    as.vector() %>% 
    data.frame() %>% 
    mutate(prediccion = case_when(
      . == 0 ~ "GTEX_B", 
      . == 1 ~ "TCGA_BLGG",
      . == 2 ~ "TCGA_GM") %>% factor()) %>% 
    select(prediccion)
  
  y_Test <- modelo %>% 
    predict(as.matrix(features_Test)) %>% 
    k_argmax() %>% 
    as.vector() %>% 
    data.frame() %>% 
    mutate(prediccion = case_when(
      . == 0 ~ "GTEX_B", 
      . == 1 ~ "TCGA_BLGG",
      . == 2 ~ "TCGA_GM") %>% factor()) %>% 
    select(prediccion)
  
  MC.Train <- table(y_Train$prediccion, Train$Clase)
  MC.Test <- table(y_Test$prediccion, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

NN3 <- function(Train, Test) {
  
  Train_Balanceado <- SMOTE(X = select(Train, -Clase), 
                          target = Train$Clase, 
                          K = 5,
                          dup_size = 2)$data %>% 
    mutate(Clase = factor(class)) %>% 
    select(-class)
  
  features_Train <- Train_Balanceado %>% 
    select(-Clase) %>% 
    data.frame()
  
  features_Test <- Test %>% 
    select(-Clase) %>% 
    data.frame()
  
  labels_Train <- Train_Balanceado %>% 
    select(Clase) %>% 
    data.frame()
  
  labels_onehot <- to_categorical(y = as.numeric(unlist(labels_Train))-1)
  
  modelo <- keras_model_sequential() %>%
    layer_dense(units = 20, 
                activation = "relu", 
                input_shape = dim(features_Train)[2]) %>% 
    layer_dense(units = 20,
                activation = "relu") %>% 
    layer_dense(units = 3, 
                activation = "softmax")
  
  modelo %>% compile(loss = "categorical_crossentropy",
                     optimizer = "adam",
                     metrics = c("accuracy"))
  
  historia <- modelo %>% 
    fit(x = as.matrix(features_Train), 
        y = labels_onehot,
        epochs = 100,
        batch_size = 32)
  
  y_Train <- modelo %>% 
    predict(as.matrix(features_Train)) %>% 
    k_argmax() %>% 
    as.vector() %>% 
    data.frame() %>% 
    mutate(prediccion = case_when(
      . == 0 ~ "GTEX_B", 
      . == 1 ~ "TCGA_BLGG",
      . == 2 ~ "TCGA_GM") %>% factor()) %>% 
    select(prediccion)
  
  y_Test <- modelo %>% 
    predict(as.matrix(features_Test)) %>% 
    k_argmax() %>% 
    as.vector() %>% 
    data.frame() %>% 
    mutate(prediccion = case_when(
      . == 0 ~ "GTEX_B", 
      . == 1 ~ "TCGA_BLGG",
      . == 2 ~ "TCGA_GM") %>% factor()) %>% 
    select(prediccion)
  
  MC.Train <- table(y_Train$prediccion, Train_Balanceado$Clase)
  MC.Test <- table(y_Test$prediccion, Test$Clase)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


M1.NN <- GenerarResultadosParalelo(Metodo = "NN1", 
                                   workers = availableCores())

M3.NN <- GenerarResultadosParalelo(Metodo = "NN3", 
                                   workers = availableCores())

