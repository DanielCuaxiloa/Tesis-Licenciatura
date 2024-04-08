
##############################################
# Modelos Gráficos Gaussianos (no dirigidos) #
# Datos Originales                           #
##############################################

rm(list = ls(all.names = TRUE))
gc()

library(bootnet)
library(qgraph)
library(readr)
library(dplyr)
library(corrplot)

## Datos en escala log2(norm_count+1)
datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

# Grupo 1, GTEX Brain (Cortex, Ba24, Ba9) ---------------------------------

DatosModelo <- datos %>% select(-c("sample","samples","_gender")) %>%
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category %in% c("GTEX Brain")) %>% 
  filter(detailed_category %in% c("Brain - Cortex",
                                  "Brain - Anterior Cingulate Cortex (Ba24)",
                                  "Brain - Frontal Cortex (Ba9)")) %>%
  select(-c("TCGA_GTEX_main_category","detailed_category"))

corMat <- cov(DatosModelo)
round(corMat, 3)

ModeloGrupo1 <- ggmModSelect(S = corMat, 
                             n = nrow(DatosModelo),
                             gamma = 0,
                             start = "full")

round(ModeloGrupo1$graph,3)

qgraph(ModeloGrupo1$graph)

ggmFit(qgraph(ModeloGrupo1$graph), corMat, nrow(DatosModelo))


Modelo <- gRim::cmod(formula = ~ .^.,
                     data = DatosModelo)

Modelo.bic <- stepwise(object = Modelo,
                       direction = "backward",
                       k = log(nrow(DatosModelo)))

round(Modelo.bic[["datainfo"]][["S"]],3)
