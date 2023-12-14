
##############################################
# Modelos Gráficos Gaussianos (no dirigidos) #
# Datos Originales                           #
##############################################

rm(list = ls(all.names = TRUE))
gc()

library(bootnet)
library(qgraph)
library(readr)
library(tidyverse)
library(corrplot)

## Datos en escala log2(norm_count+1)
Datos <- read_tsv("../Datos/denseDataOnlyDownload.tsv")

Datos <- Datos %>% select(-c("sample","samples")) %>%
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category %in% c("GTEX Brain", 
                                        "TCGA Brain Lower Grade Glioma",
                                        "TCGA Glioblastoma Multiforme")) %>% 
  filter(detailed_category %in% c("Brain - Cortex",
                                  "Brain - Anterior Cingulate Cortex (Ba24)",
                                  "Brain - Frontal Cortex (Ba9)",
                                  "Brain Lower Grade Glioma",
                                  "Glioblastoma Multiforme")) %>%
  mutate(TCGA_GTEX_main_category = fct_relevel(droplevels(TCGA_GTEX_main_category),
                                               c("GTEX Brain", 
                                                 "TCGA Brain Lower Grade Glioma",
                                                 "TCGA Glioblastoma Multiforme"))) %>% 
  select(-c("_gender", "detailed_category"))

Datos <- Datos %>% 
  rename(y = TCGA_GTEX_main_category) %>% 
  mutate(Clase = factor(case_when(y == "GTEX Brain" ~ "GTEX_B",
                                  y == "TCGA Brain Lower Grade Glioma" ~ "TCGA_BLGG",
                                  y == "TCGA Glioblastoma Multiforme" ~ "TCGA_GM"))) %>% 
  select(-y)

# Grupo 1, GTEX Brain (Cortex, Ba24, Ba9) ---------------------------------

Grupo1 <- Datos %>% 
  filter(Clase == "GTEX_B")

Grupo1 <- model.matrix(object = Clase~.^2,
                       data = Grupo1)

Grupo1 <- data.frame(Grupo1) %>% 
  select(-c("X.Intercept."))

ModeloGrupo1 <- estimateNetwork(data = Grupo1,
                                default = "EBICglasso",
                                corMethod = "cor",
                                tuning = 5)
plot(ModeloGrupo1,
     edge.labels = FALSE,
     font = 2)


# Grupo 2, TCGA Brain Lower Grade Glioma ----------------------------------

Grupo2 <- Datos %>% 
  filter(Clase == "TCGA_BLGG")

Grupo2 <- model.matrix(object = Clase~.^2,
                       data = Grupo2)

Grupo2 <- data.frame(Grupo2) %>% 
  select(-c("X.Intercept."))

ModeloGrupo2 <- estimateNetwork(data = Grupo2,
                                default = "EBICglasso",
                                corMethod = "cor",
                                tuning = 5)
plot(ModeloGrupo2,
     edge.labels = FALSE,
     font = 2)


# Grupo 3, TCGA Glioblastoma Multiforme -----------------------------------

Grupo3 <- Datos %>% 
  filter(Clase == "TCGA_GM")

Grupo3 <- model.matrix(object = Clase~.^2,
                       data = Grupo3)

Grupo3 <- data.frame(Grupo3) %>% 
  select(-c("X.Intercept."))

ModeloGrupo3 <- estimateNetwork(data = Grupo3,
                                default = "EBICglasso",
                                corMethod = "cor",
                                tuning = 2.5)
plot(ModeloGrupo3,
     edge.labels = FALSE,
     font = 2)



