
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

Modelo6 <- read.csv(file = "Modelo6.csv", 
                     stringsAsFactors = TRUE)


Datos <- bind_rows(Modelo1, Modelo2, Modelo3, Modelo4, Modelo5, Modelo6) %>% 
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
  arrange(desc(TestClase2)) %>% 
  ungroup()


