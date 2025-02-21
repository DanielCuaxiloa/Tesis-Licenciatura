
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
            size = 4) +
scale_fill_gradient(low = "white", high = "steelblue1") + 
facet_grid(cols = vars(Modelo),
             rows = vars(Nombre)) + 
  labs(title = "Matrices de Confusión",
       subtitle = "Promedios",
       x = "Predicción", 
       y = "Referencia") +
  theme_bw() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

