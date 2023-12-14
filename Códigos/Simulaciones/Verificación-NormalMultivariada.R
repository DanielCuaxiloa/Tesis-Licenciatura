
################################
# Verificación del supuesto de # 
# Normalidad Multivariada      #      
################################

rm(list = ls(all.names = TRUE))
gc()

library(readr)
library(dplyr)
library(tidyr)

library(gtools)

library(MVN)
library(mvnormtest)
library(mnt)
library(energy)

library(huge)

library(bootnet)

## Datos en escala log2(norm_count+1)
datos <- read_tsv("../Datos/denseDataOnlyDownload.tsv")


# Grupo 1, GTEX Brain (Cortex, Ba24, Ba9) ---------------------------------

Grupo1 <- datos %>% select(-c("sample","samples","_gender")) %>%
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category %in% c("GTEX Brain")) %>% 
  filter(detailed_category %in% c("Brain - Cortex",
                                  "Brain - Anterior Cingulate Cortex (Ba24)",
                                  "Brain - Frontal Cortex (Ba9)")) %>%
  select(-c("TCGA_GTEX_main_category","detailed_category"))

## Estadísticas descriptivas
Grupo1 %>% pivot_longer(cols = everything(),
                        names_to = "GeneExpression",
                        values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo1, mvnTest = "royston")
# mshapiro.test(U = as.matrix(Grupo1))

# TRUE -> Evidencia en contra.
test.BHEP(data = as.matrix(Grupo1), 
          MC.rep = 100, 
          alpha = 0.05)[["Decision"]]

mvnorm.test(x = as.matrix(Grupo1), 
            R = 100)

## Transformación
Grupo1T <- huge.npn(Grupo1, 
                    npn.func = "shrinkage")

## Estadísticas descriptivas
Grupo1T %>% data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = "GeneExpression",
               values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo1T, mvnTest = "royston")
mvn(data = Grupo1T, mvnTest = "hz")
mvn(data = Grupo1T, mvnTest = "royston")
# mshapiro.test(Grupo1T)
test.BHEP(data = Grupo1T, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]
mvnorm.test(x = Grupo1T, 
            R = 1000)


# Grupo 2, TCGA Brain Lower Grade Glioma ----------------------------------

Grupo2 <- datos %>% select(-c("sample","samples","_gender")) %>%
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category %in% c("TCGA Brain Lower Grade Glioma")) %>% 
  select(-c("TCGA_GTEX_main_category","detailed_category"))

## Estadísticas descriptivas
Grupo2 %>% pivot_longer(cols = everything(),
                        names_to = "GeneExpression",
                        values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo2, mvnTest = "royston")
# mshapiro.test(U = as.matrix(Grupo2))
test.BHEP(data = as.matrix(Grupo2), 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]
mvnorm.test(x = as.matrix(Grupo2), 
            R = 100)

## Transformación
Grupo2T <- huge.npn(Grupo2, 
                    npn.func = "shrinkage")

## Estadísticas descriptivas
Grupo2T %>% data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = "GeneExpression",
               values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo2T, mvnTest = "mardia")
mvn(data = Grupo2T, mvnTest = "hz")
mvn(data = Grupo2T, mvnTest = "royston")
# mshapiro.test(Grupo2T)
test.BHEP(data = Grupo2T, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]
mvnorm.test(x = Grupo2T, 
            R = 100)


# Grupo 3, TCGA Glioblastoma Multiforme -----------------------------------

Grupo3 <- datos %>% select(-c("sample","samples","_gender")) %>%
  mutate(TCGA_GTEX_main_category = factor(TCGA_GTEX_main_category),
         detailed_category = factor(detailed_category)) %>% 
  filter(TCGA_GTEX_main_category %in% c("TCGA Glioblastoma Multiforme")) %>%
  select(-c("TCGA_GTEX_main_category","detailed_category"))

## Estadísticas descriptivas
Grupo3 %>% pivot_longer(cols = everything(),
                        names_to = "GeneExpression",
                        values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo3, mvnTest = "royston")
# mshapiro.test(U = as.matrix(Grupo3))
test.BHEP(data = as.matrix(Grupo3), 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]
mvnorm.test(x = as.matrix(Grupo3), 
            R = 100)

## Transformación
Grupo3T <- huge.npn(Grupo3, npn.func = "shrinkage")

## Estadísticas descriptivas
Grupo3T %>% data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = "GeneExpression",
               values_to = "Aux") %>% 
  group_by(GeneExpression) %>% 
  summarise(Minimo = min(Aux),
            Maximo = max(Aux),
            Media = mean(Aux),
            Varianza = var(Aux)) %>% 
  arrange(Media)

mvn(data = Grupo3T, mvnTest = "mardia")
mvn(data = Grupo3T, mvnTest = "hz")
mvn(data = Grupo3T, mvnTest = "royston")
# mshapiro.test(Grupo3T)
test.BHEP(data = Grupo3T, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]
mvnorm.test(x = Grupo3T, 
            R = 1000)

# Algoritmo 1 -------------------------------------------------------------

v <- seq(1:11)
r <- 11
B <- FALSE
while (B == FALSE) {
  
  indices <- combinations(n = 11, r = r, v = v)
  
  for (i in 1:dim(indices)[1]) {
    
    Grupo1Aux <- Grupo1T[,indices[i,]]
    Grupo2Aux <- Grupo2T[,indices[i,]]
    Grupo3Aux <- Grupo3T[,indices[i,]]
    
    Prueba_Grupo1Aux <- mvn(data = Grupo1Aux, mvnTest = "royston")
    Prueba_Grupo2Aux <- mvn(data = Grupo2Aux, mvnTest = "royston")
    Prueba_Grupo3Aux <- mvn(data = Grupo3Aux, mvnTest = "royston")
    
    if(Prueba_Grupo1Aux[["multivariateNormality"]][["MVN"]] == "YES" &
       Prueba_Grupo2Aux[["multivariateNormality"]][["MVN"]] == "YES" &
       Prueba_Grupo3Aux[["multivariateNormality"]][["MVN"]] == "YES") {
      
      B <- TRUE
      Variables <- colnames(Grupo1Aux)
      break
      
    } else {
      
      B <- FALSE
      
    }
    
  }       
  
  r <- r-1
  
}

mvn(data = Grupo1T[,Variables], mvnTest = "royston")$multivariateNormality
mvn(data = Grupo2T[,Variables], mvnTest = "royston")$multivariateNormality
mvn(data = Grupo3T[,Variables], mvnTest = "royston")$multivariateNormality

