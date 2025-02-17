
##############################
# Recopilación de Resultados #
##############################

library(dplyr)
library(ggplot2)

Modelo1 <- read.csv(file = "Modelo1.csv", 
                     stringsAsFactors = TRUE)

Modelo2 <- read.csv(file = "Modelo2.csv", 
                     stringsAsFactors = TRUE)

Modelo3 <- read.csv(file = "Modelo3.csv", 
                     stringsAsFactors = TRUE)

Modelo4 <- read.csv(file = "Modelo4.csv", 
                     stringsAsFactors = TRUE)

Modelo5 <- read.csv(file = "Modelo5.csv", 
                     stringsAsFactors = TRUE)

Datos <- bind_rows(Modelo1, Modelo2, Modelo3, Modelo4, Modelo5) %>% 
  select(-ID1) %>% 
  select(Modelo, Nombre, TestClase1, TestClase2, TestClase3, TestGlobal) %>% 
  rename(
    'TCC General' = TestGlobal,
    'TCC - GTEX Brain' = TestClase1,
    'TCC - TCGA LGG' = TestClase2,
    'TCC - TCGA GM' = TestClase3,
    Model = Nombre
  )
  

ggplot(data = Datos,
       mapping = aes(x = Model, y = `TCC - TCGA GM`)) +
  facet_grid(cols = vars(Modelo), 
             scales = "free_x") +
  geom_boxplot(fill = "steelblue3") + 
  labs(x = "Modelo", 
       title = "Comparación de modelos") +
  theme_bw()

Datos %>% 
  group_by(Modelo, Model) %>% 
  summarise(TCC_G = mean(`TCC General`),
            TCC_C1 = mean(`TCC - GTEX Brain`),
            TCC_C2 = mean(`TCC - TCGA LGG`),
            TCC_C3 = mean(`TCC - TCGA GM`)) %>% 
  arrange(desc(TCC_C1)) %>% 
  ungroup()


