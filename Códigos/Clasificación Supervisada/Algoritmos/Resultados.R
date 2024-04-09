
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

Datos <- bind_rows(Modelo_1, Modelo_2, Modelo_3) %>% 
  select(-ID1)

ggplot(data = Datos,
       mapping = aes(x = Nombre, y = TestGlobal)) +
  facet_wrap(~Modelo, scales = "free_x") +
  geom_boxplot(fill = "steelblue3") + 
  theme_bw()

