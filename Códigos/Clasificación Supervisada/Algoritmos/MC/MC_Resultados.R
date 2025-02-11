
#################################
# Recopilación de Resultados MC #
#################################

library(dplyr)
library(forcats)
library(ggplot2)

Modelo1 <- read.csv(file = "MC1.csv", 
                    stringsAsFactors = TRUE)

Modelo2 <- read.csv(file = "MC2.csv", 
                    stringsAsFactors = TRUE)

Modelo3 <- read.csv(file = "MC3.csv", 
                    stringsAsFactors = TRUE)

Modelo4 <- read.csv(file = "MC4.csv", 
                    stringsAsFactors = TRUE)

Modelo5 <- read.csv(file = "MC5.csv", 
                    stringsAsFactors = TRUE)

Datos <- bind_rows(Modelo1, Modelo2, Modelo3, Modelo4, Modelo5) %>%
  mutate(Referencia = fct_relevel(droplevels(Referencia),
                                  c("TCGA GM",
                                    "TCGA LGG",
                                    "GTEX Brain")),
         Prediccion = fct_relevel(droplevels(Prediccion),
                                  c("GTEX Brain",
                                    "TCGA LGG",
                                    "TCGA GM")))

ggplot(Datos, 
       aes(x = Prediccion, y = Referencia, fill = Porcentaje)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f%%", Porcentaje)), 
            color = "black", 
            size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_grid(cols = vars(Modelo),
             rows = vars(Nombre)) + 
  labs(title = "Promedio Matrices de Confusión", x = "Predicción", y = "Referencia") +
  theme_bw()
