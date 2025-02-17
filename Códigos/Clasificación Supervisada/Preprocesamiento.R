
###############################
# Preprocesamiento | 3 Clases # 
###############################

rm(list = ls(all.names = TRUE))
gc()

library(readr)
library(dplyr)
library(forcats)
library(rsample)

# Datos -------------------------------------------------------------------

## Lectura de datos originales.
Datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

## Pre-procesamiento de datos originales.
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
  select(-c("_gender", "detailed_category")) %>% 
  rename(y = TCGA_GTEX_main_category) %>% 
  mutate(Clase = factor(case_when(y == "GTEX Brain" ~ "GTEX_B",
                                  y == "TCGA Brain Lower Grade Glioma" ~ "TCGA_BLGG",
                                  y == "TCGA Glioblastoma Multiforme" ~ "TCGA_GM"))) %>% 
  select(-y)

summary(Datos)

## Remuestreo de datos en conjuntos de entrenamiento (Train) y evaluacion (Test)
## V-K Cross-validation (V = 20, K = 5).
set.seed(1234)
Folds <- vfold_cv(data = Datos,
                  v = 5,
                  repeats = 100,
                  strata = Clase)

rm(Datos)

save.image("Folds.RData")






