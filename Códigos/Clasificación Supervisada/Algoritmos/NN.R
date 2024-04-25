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
  
  MC.Train <- table(Train$Clase, y_Train$prediccion)
  MC.Test <- table(Test$Clase, y_Test$prediccion)
  
  return(list(MC.Train = MC.Train,
              MC.Test = MC.Test))
  
}

M1.NN <- Evaluacion(Metodo = "NN1",
                    workers = availableCores())

# Gráficas ----------------------------------------------------------------

G.NN_1 <- M1.NN[["Global"]] %>% 
  mutate(Modelo = "M5",
         Nombre = "NN_1")

M5 <- bind_rows(G.NN_1) %>% 
  mutate(Modelo = as.factor(Modelo),
         Nombre = as.factor(Nombre))

ggplot(data = M5,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

mean(G.NN_1$TestGlobal)

write.csv(x = M5,
          file = "Modelo_5.csv",
          row.names = FALSE)



