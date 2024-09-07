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

library(keras)
library(tensorflow)
library(kerastuneR)

library(ggplot2)

# Funciones Auxiliares ----------------------------------------------------

source("../FuncionesAuxiliares.R")

# Conjuntos Train y Test --------------------------------------------------

load("../Folds.RData")


# Esquemas de clasificación -----------------------------------------------

set.seed(123)
set_random_seed(123)

NN.1 <- function(Train, Test) {
  
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
  
  MC.Train <- table(Train$Clase, y_Train$prediccion)
  MC.Test <- table(Test$Clase, y_Test$prediccion)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

NN.2 <- function(Train, Test) {
  
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
  
  MC.Train <- table(Train$Clase, y_Train$prediccion)
  MC.Test <- table(Test$Clase, y_Test$prediccion)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}


# Resultados --------------------------------------------------------------

M6.NN.1 <- Evaluacion(Metodo = "NN.1",
                    workers = availableCores())

M6.NN.2 <- Evaluacion(Metodo = "NN.2",
                    workers = availableCores())

# Gráficas ----------------------------------------------------------------

G.NN.1 <- M6.NN.1[["Global"]] %>% 
  mutate(Modelo = "M6",
         Nombre = "NN.1")

G.NN.2 <- M6.NN.2[["Global"]] %>% 
  mutate(Modelo = "M6",
         Nombre = "NN.2")

M6 <- bind_rows(G.NN.1, G.NN.2) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M6,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.NN.1$TestGlobal)
mean(G.NN.2$TestGlobal)


write.csv(x = M6,
          file = "Modelo6.csv",
          row.names = FALSE)



