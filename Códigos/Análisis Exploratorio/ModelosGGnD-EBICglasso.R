
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
Datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

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

GTEX_B <- Datos %>% 
  filter(Clase == "GTEX_B") %>% 
  select(-Clase)

## Correlación de Pearson
ModeloGrupo1.1 <- estimateNetwork(data = GTEX_B,
                                  default = "cor")

## EBICglasso
ModeloGrupo1.2 <- estimateNetwork(data = GTEX_B,
                                  default = "EBICglasso",
                                  corMethod = "cor",
                                  tuning = 0.5)

plot(ModeloGrupo1.1,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

plot(ModeloGrupo1.2,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

# Grupo 2, TCGA Brain Lower Grade Glioma ----------------------------------

TCGA_BLGG <- Datos %>% 
  filter(Clase == "TCGA_BLGG") %>% 
  select(-Clase)

## Correlción de Pearson
ModeloGrupo2.1 <- estimateNetwork(data = TCGA_BLGG,
                                  default = "cor")

## EBICglasso
ModeloGrupo2.2 <- estimateNetwork(data = TCGA_BLGG,
                                  default = "EBICglasso",
                                  corMethod = "cor",
                                  tuning = 0.5)

plot(ModeloGrupo2.1,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

plot(ModeloGrupo2.2,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

# Grupo 3, TCGA Glioblastoma Multiforme -----------------------------------

TCGA_GM <- Datos %>% 
  filter(Clase == "TCGA_GM") %>% 
  select(-Clase)

ModeloGrupo3.1 <- estimateNetwork(data = TCGA_GM,
                                 default = "cor")

ModeloGrupo3.2 <- estimateNetwork(data = TCGA_GM,
                                  default = "EBICglasso",
                                  corMethod = "cor",
                                  tuning = 0.5)

plot(ModeloGrupo3.1,
     layout = "circle",
     edge.labels = FALSE,
     font = 2) 

plot(ModeloGrupo3.2,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)
