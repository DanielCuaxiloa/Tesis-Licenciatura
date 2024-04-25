
##############################
# Recopilación de Resultados #
##############################

library(dplyr)
library(ggplot2)

Modelo_1 <- read.csv(file = "Modelo_1.csv", 
                     stringsAsFactors = TRUE)

Modelo_2 <- read.csv(file = "Modelo_2.csv", 
                     stringsAsFactors = TRUE)

Modelo_3 <- read.csv(file = "Modelo_3.csv", 
                     stringsAsFactors = TRUE)

Modelo_4 <- read.csv(file = "Modelo_4.csv", 
                     stringsAsFactors = TRUE)

Modelo_5 <- read.csv(file = "Modelo_5.csv", 
                     stringsAsFactors = TRUE)

Datos <- bind_rows(Modelo_1, Modelo_2, Modelo_3, Modelo_4, Modelo_5) %>% 
  select(-ID1)

ggplot(data = Datos,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  facet_grid(cols = vars(Modelo), 
             scales = "free_x") +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

Datos %>% 
  group_by(Modelo, Nombre) %>% 
  summarise(TestGlobal = mean(TestGlobal),
            TestClase1 = mean(TestClase1),
            TestClase2 = mean(TestClase2),
            TestClase3 = mean(TestClase3)) %>% 
  arrange(desc(TestGlobal)) %>% 
  ungroup()


