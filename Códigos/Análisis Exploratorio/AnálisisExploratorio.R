
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
library(ggplot2)

library(bootnet)

# Datos -------------------------------------------------------------------

## Datos en escala log2(norm_count+1)
Datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

## Preprocesamiento de datos originales.
Datos <- Datos %>% 
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
  select(-c("_gender", "detailed_category")) %>% 
  mutate(Clase = factor(case_when(TCGA_GTEX_main_category == "GTEX Brain" ~ "GTEX_Brain",
                                  TCGA_GTEX_main_category == "TCGA Brain Lower Grade Glioma" ~ "TCGA_LGG",
                                  TCGA_GTEX_main_category == "TCGA Glioblastoma Multiforme" ~ "TCGA_GM"))) %>% 
  select(-TCGA_GTEX_main_category)

write.csv(Datos, file = "Datos.csv", row.names = FALSE)

# Distancias --------------------------------------------------------------

D <- Datos %>% 
  pivot_longer(cols = -Clase,
               names_to = "GeneExpression",
               values_to = "Aux") %>% 
  group_by(Clase, GeneExpression) %>% 
  summarise(Media = mean(Aux)) %>% 
  ungroup() %>%
  pivot_wider(names_from = Clase, values_from = Media) %>%
  mutate(across(starts_with("TCGA_"), ~ abs(. - GTEX_Brain))) %>% 
  select(GeneExpression, starts_with("TCGA_")) %>%
  pivot_longer(cols = starts_with("TCGA_"), names_to = "Clase", values_to = "Diferencia")

D$Clase <- factor(D$Clase, levels=c("GTEX_Brain", "TCGA_LGG", "TCGA_GM"))

  ggplot(data = D,
        mapping = aes(x = Clase, y = reorder(x = GeneExpression,
                                    X = Diferencia,
                                    FUN = median), fill = Diferencia)) +
  geom_tile() +
  geom_text(aes(label = round(Diferencia, 1)), color = "steelblue1") +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Grupo",
       y = "Gene Expression",
       fill = "")

# Gráfica BoxPlot por grupos ----------------------------------------------

Datos %>% 
pivot_longer(cols = -Clase,
             names_to = "GeneExpression",
             values_to = "Aux") %>% 
ggplot(mapping = aes(x = Aux, 
                     y = reorder(x = GeneExpression,
                                 X = Aux,
                                 FUN = median),
                     fill = Clase)) + 
  geom_boxplot(alpha = 0.5) +
  scale_fill_viridis_d() +
  facet_wrap(~Clase) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(x = "log2(norm count + 1)",
       y = "Enzima")


# Gráfica de densidades por grupos ----------------------------------------

Datos %>% 
  pivot_longer(cols = -Clase,
               names_to = "GeneExpression",
               values_to = "Aux") %>% 
  ggplot(mapping = aes(x = Aux, 
                       y = reorder(x = GeneExpression,
                                   X = Aux,
                                   FUN = median),
                       fill = Clase)) + 
  geom_density_ridges(alpha = 0.5,
                      scale = 0.9) +
  scale_fill_viridis_d() +
  facet_wrap(~Clase) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(x = "log2(norm count + 1)",
       y = "Enzima")


# Correlación de Pearson --------------------------------------------------

## Clase: Brain_Cortex: * Brain - Cortex
##                      * Brain - Anterior Cingulate Cortex (Ba24)
##                      * Brain - Frontal Cortex (Ba9))
Pearson.GTEX_Brain <- Datos %>% 
  filter(Clase == "GTEX_Brain") %>%
  select(-Clase) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  data.frame()

MGG.Pearson.GTEX_Brain <- Datos %>% 
  filter(Clase == "GTEX_Brain") %>%
  select(-Clase) %>% 
  estimateNetwork(default = "cor")

## Clase: TCGA Lower Grade Glioma
Pearson.TCGA_LGG <- Datos %>% 
  filter(Clase == "TCGA_LGG") %>%
  select(-Clase) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  data.frame()

MGG.Pearson.TCGA_LGG <- Datos %>% 
  filter(Clase == "TCGA_LGG") %>%
  select(-Clase) %>% 
  estimateNetwork(default = "cor")

## Clase: TCGA Glioblastoma Multiforme
Pearson.TCGA_GBM <- Datos %>%
  filter(Clase == "TCGA_GM") %>%
  select(-Clase) %>% 
  cor(method = "pearson") %>% 
  round(digits = 3) %>% 
  data.frame()

MGG.Pearson.TCGA_GM <- Datos %>% 
  filter(Clase == "TCGA_GM") %>%
  select(-Clase) %>% 
  estimateNetwork(default = "cor")


# GGM Gráficas ------------------------------------------------------------

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
plot(MGG.Pearson.GTEX_Brain,
     title = "GTEX_Brain",
     layout = "circle",
     edge.labels = FALSE,
     label.prop = 1,
     font = 2)
plot(MGG.Pearson.TCGA_LGG,
     title = "TCGA_LGG",
     layout = "circle",
     edge.labels = FALSE,
     label.prop = 1,
     font = 2)
plot(MGG.Pearson.TCGA_GM,
     title = "TCGA_GM",
     layout = "circle",
     edge.labels = FALSE,
     label.prop = 1,
     font = 2)




