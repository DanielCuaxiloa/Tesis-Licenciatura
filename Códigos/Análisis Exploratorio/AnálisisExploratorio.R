
#########################
# Análisis Exploratorio #
#########################

rm(list = ls(all.names = TRUE))
gc()

library(readr)
library(tidyverse)
library(ggforce)
library(ggridges)
library(corrplot)


# Datos -------------------------------------------------------------------

## Datos en escala log2(norm_count+1)
datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

## Pre-procesamiento de datos originales para gráfica BoxPlot por grupos
datos_aux <- datos %>% 
  select(-c("sample","samples")) %>%
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
  pivot_longer(cols = !c("_gender", "TCGA_GTEX_main_category", "detailed_category"),
               names_to = "GeneExpression",
               values_to = "Aux")


# Gráfica BoxPlot por grupos ----------------------------------------------

ggplot(data = datos_aux,
       mapping = aes(x = reorder(x = GeneExpression,
                                 X = Aux,
                                 FUN = median), 
                     y = Aux,
                     fill = TCGA_GTEX_main_category)) + 
  geom_boxplot(alpha = 0.5) +
  scale_fill_viridis_d() +
  facet_wrap(~TCGA_GTEX_main_category) +
  coord_flip() + 
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8)) + 
  labs(x="Gene Expression",
       y="log2(norm count + 1)")


# Gráfica de densidades por grupos ----------------------------------------

ggplot(data = datos_aux,
       mapping = aes(y = reorder(x = GeneExpression,
                                 X = Aux,
                                 FUN = median), 
                     x = Aux,
                     fill = TCGA_GTEX_main_category)) + 
  geom_density_ridges(alpha = 0.5,
                      scale = 0.9) +
  scale_fill_viridis_d() +
  facet_wrap(~TCGA_GTEX_main_category) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8)) + 
  labs(y="Gene Expression",
       x="log2(norm count + 1)")


# Correlación de Pearson --------------------------------------------------

## Grupo: GTEX Brain(Brain - Cortex, Brain - Anterior Cingulate Cortex (Ba24) y
## Brain - Frontal Cortex (Ba9))
datos %>% 
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category == "GTEX Brain") %>% 
  filter(detailed_category %in% c("Brain - Cortex",
                                  "Brain - Anterior Cingulate Cortex (Ba24)",
                                  "Brain - Frontal Cortex (Ba9)")) %>%
  select(!c("sample","samples","_gender",
            "detailed_category","TCGA_GTEX_main_category")) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  corrplot(method = "ellipse",
           type = "upper",
           title = "GTEX Brain \n correlación de Pearson",
           mar=c(0,0,2,0), 
           diag = FALSE,
           addCoef.col = "black",
           tl.cex = 0.75)

## Grupo: TCGA Brain Lower Grade Glioma
datos %>% 
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category)) %>% 
  filter(TCGA_GTEX_main_category == "TCGA Brain Lower Grade Glioma") %>%
  select(!c("sample","samples","_gender",
            "detailed_category","TCGA_GTEX_main_category")) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  corrplot(method = "ellipse",
           type = "upper",
           title = "TCGA Brain Lower Grade Glioma \n correlación de Pearson",
           mar=c(0,0,2,0), 
           diag = FALSE,
           addCoef.col = "black",
           tl.cex = 0.75)

## Grupo: TCGA Glioblastoma Multiforme
datos %>% 
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category)) %>% 
  filter(TCGA_GTEX_main_category == "TCGA Glioblastoma Multiforme") %>% 
  select(!c("sample","samples","_gender",
            "detailed_category","TCGA_GTEX_main_category")) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  corrplot(method = "ellipse",
           type = "upper",
           title = "TCGA Glioblastoma Multiforme \n Correlación de Pearson",
           mar=c(0,0,2,0), 
           diag = FALSE,
           addCoef.col = "black",
           tl.cex = 0.75)






