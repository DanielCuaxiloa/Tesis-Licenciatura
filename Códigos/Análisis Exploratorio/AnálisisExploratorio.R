
#########################
# Análisis Exploratorio #
#########################


rm(list = ls(all.names = TRUE))
gc()

library(readr)
library(ppcor)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggcorrplot)
library(rstatix)
library(ggpubr)
library(GGally)
library(cowplot)
library(gridExtra)
library(ggbiplot)


# Datos -------------------------------------------------------------------

## Datos en escala log2(norm_count+1)
Datos <- read_tsv("../../Datos/denseDataOnlyDownload.tsv")

## Preprocesamiento de datos originales.
Datos <- Datos %>% 
  dplyr::select(-c("sample","samples")) %>%
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
  dplyr::select(-c("_gender", "detailed_category")) %>% 
  mutate(Clase = factor(case_when(TCGA_GTEX_main_category == "GTEX Brain" ~ "GTEX_Brain",
                                  TCGA_GTEX_main_category == "TCGA Brain Lower Grade Glioma" ~ "TCGA_LGG",
                                  TCGA_GTEX_main_category == "TCGA Glioblastoma Multiforme" ~ "TCGA_GM"))) %>% 
  dplyr::select(-TCGA_GTEX_main_category)

## write.csv(Datos, file = "Datos.csv", row.names = FALSE)


# Correlación de Pearson --------------------------------------------------

## Clase: Brain_Cortex: * Brain - Cortex
##                      * Brain - Anterior Cingulate Cortex (Ba24)
##                      * Brain - Frontal Cortex (Ba9))
Pearson.GTEX_Brain <- Datos %>% 
  filter(Clase == "GTEX_Brain") %>%
  dplyr::select(-Clase) %>% 
  pcor(method = "pearson")

Pearson.GTEX_Brain <- Pearson.GTEX_Brain$estimate %>%
  round(digits = 2) %>% 
  data.frame()

ggcorrplot(Pearson.GTEX_Brain, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           type = "lower", 
           lab = TRUE, 
           lab_size = 5,
           lab_col = "black") +
  labs(title = "Matriz de correlaciones parciales",
       subtitle = "Grupo GTEX Brain",
       x = "",
       y = "") + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

## Clase: TCGA Lower Grade Glioma
Pearson.TCGA_LGG <- Datos %>% 
  filter(Clase == "TCGA_LGG") %>%
  dplyr::select(-Clase) %>% 
  pcor(method = "pearson")

Pearson.TCGA_LGG <- Pearson.TCGA_LGG$estimate %>%
  round(digits = 2) %>% 
  data.frame()

ggcorrplot(Pearson.TCGA_LGG, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           type = "lower", 
           lab = TRUE, 
           lab_size = 5,
           lab_col = "black") +
  labs(title = "Matriz de correlaciones parciales",
       subtitle = "Grupo TCGA LGG",
       x = "",
       y = "") + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

## Clase: TCGA Glioblastoma Multiforme
Pearson.TCGA_GM <- Datos %>%
  filter(Clase == "TCGA_GM") %>%
  dplyr::select(-Clase) %>% 
  pcor(method = "pearson")

Pearson.TCGA_GM <- Pearson.TCGA_GM$estimate %>% 
  round(digits = 2) %>% 
  data.frame()

ggcorrplot(Pearson.TCGA_GM, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           type = "lower", 
           lab = TRUE, 
           lab_size = 5,
           lab_col = "black") +
  labs(title = "Matriz de correlaciones parciales",
       subtitle = "Grupo TCGA GM",
       x = "",
       y = "") + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


# Diagrama de dispersión --------------------------------------------------

Datos <- Datos %>% 
  mutate(Clase = factor(recode(Clase,
                               "GTEX_Brain" = "GTEX Brain",
                               "TCGA_LGG"   = "TCGA LGG",
                               "TCGA_GM"    = "TCGA GM"),
                        levels = c("GTEX Brain", "TCGA LGG", "TCGA GM")))

Aux1 <- ggpairs(Datos,
                mapping = aes(color = Clase),
                columns = 1:11,
                lower = list(continuous = wrap("points", alpha = 0.7, size = 0.5)), 
                #diag  = list(continuous = wrap("cor", size = 5, mapping = aes())), 
                upper = list(continuous = wrap("cor", size = 3, 
                                               method = "pearson", fontface="bold")),
                legend = 1) +
  scale_color_brewer(palette = "Set2", name = "Clase") +
  labs(title = "Diagramas de dispersión",
       subtitle = "Principales enzimas involucradas en la vía de las Kinureninas") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

Aux2 <- ggplot(Datos, aes(x = 1, y = 1, color = Clase)) +
  geom_point() +
  scale_color_brewer(palette = "Set2", name = "Clase") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

legenda <- get_legend(Aux2)
Diagrama <- ggmatrix_gtable(Aux1) 

plot_grid(Diagrama, legenda, rel_widths = c(6, 1))


# Componentes principales -------------------------------------------------

pca <- prcomp(Datos[, -which(names(Datos) == "Clase")], 
              center = TRUE, 
              scale. = FALSE)

var_exp <- summary(pca)$importance[2, 1:3] * 100

## PC1 v.s PC2
ggbiplot(pca, group = Datos$Clase, 
         choices = c(1, 2),
         ellipse = TRUE, 
         ellipse.alpha = 0.0) +
  theme_bw() +
  labs(title = "Componentes Principales", 
       subtitle = "Análisis de la varianza explicada por los componentes",
       x = paste0("PC1 (", round(var_exp[1], 2), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 2), "%)")) +
  scale_color_brewer(palette = "Set2", name = "Clase") +
  guides(fill = "none", 
         color = guide_legend(title = "Clase")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

## PC1 v.s PC3

ggbiplot(pca, group = Datos$Clase, 
         choices = c(1, 3),
         ellipse = TRUE, 
         ellipse.alpha = 0.05) +
  theme_bw() +
  labs(title = "Componentes Principales", 
       subtitle = "Análisis de la varianza explicada por los componentes",
       x = paste0("PC1 (", round(var_exp[1], 2), "%)"),
       y = paste0("PC3 (", round(var_exp[3], 2), "%)")) +
  scale_color_brewer(palette = "Set2", name = "Clase") +
  guides(fill = "none", 
         color = guide_legend(title = "Clase")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

## PC2 v.s PC3

ggbiplot(pca, group = Datos$Clase, 
         choices = c(2, 3),
         ellipse = TRUE, 
         ellipse.alpha = 0.05) +
  theme_bw() +
  labs(title = "Componentes Principales", 
       subtitle = "Análisis de la varianza explicada por los componentes",
       x = paste0("PC2 (", round(var_exp[2], 2), "%)"),
       y = paste0("PC3 (", round(var_exp[3], 2), "%)")) +
  scale_color_brewer(palette = "Set2", name = "Clase") +
  guides(fill = "none", 
         color = guide_legend(title = "Clase")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) 

PCA <- pca$x %>% 
  as.data.frame() %>% 
  select(PC1, PC2, PC3) %>% 
  bind_cols(Datos) %>% 
  select_if(is.numeric) %>% 
  cor(use = "complete.obs") %>%
  select(PC1, PC2, PC3)

round(PCA[, c("PC1", "PC2", "PC3")],2)


# Pruebas Kruskal-Wallis y Dunn ---------------------------------------

## Lista de resultados
dunn_results <- list()
kw <- data.frame()

## Comparaciones KW y Dunn
for (var in names(Datos)[1:11]) {
  
  formula <- as.formula(paste(var, "~ Clase"))
  
  # Prueba de Kruskal-Wallis
  kw_test <- kruskal_test(Datos, formula = formula)
  
  kw <- rbind(kw, data.frame(variable = var, 
                             p.value = kw_test$p))
  
  # Prueba de Dunn con corrección Bonferroni
  dunn_test <- Datos %>%
    dunn_test(formula = formula, p.adjust.method = "bonferroni")
  
  dunn_test$variable <- var
  
  # Ajustar posiciones para etiquetas de p-values
  max_y <- max(Datos[[var]], na.rm = TRUE)
  step <- 0.1 * max_y
  
  dunn_test <- dunn_test %>%
    mutate(y.position = seq(max_y + step, by = step, length.out = nrow(dunn_test)),
           group1 = recode(group1, 
                           "GTEX_Brain" = "GTEX Brain", 
                           "TCGA_LGG" = "TCGA LGG", 
                           "TCGA_GM" = "TCGA GM"),
           group2 = recode(group2, 
                           "GTEX_Brain" = "GTEX Brain", 
                           "TCGA_LGG" = "TCGA LGG", 
                           "TCGA_GM" = "TCGA GM"))
  
  dunn_results[[var]] <- dunn_test
}

dunn_df <- do.call(rbind, dunn_results)

# Transformar los datos a formato largo
Datos_Long <- Datos %>%
  pivot_longer(cols = names(Datos)[1:11], 
               names_to = "variable", 
               values_to = "value")

ggplot(Datos_Long, 
       aes(x = Clase, y = value)) +
  geom_jitter(aes(color = Clase), 
              width = 0.2, 
              alpha = 0.5) +
  geom_boxplot(aes(fill = Clase), 
               alpha = 0.6, 
               outlier.shape = NA) +
  stat_pvalue_manual(dunn_df, 
                     label = "p.adj.signif", 
                     hide.ns = FALSE, 
                     size = 3) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Clase", 
       y = "log2(norm count + 1)", 
       title = "Kruskal-Wallis Test",
       subtitle = "Dunn's Test | Post Hoc") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(data = kw, 
            aes(x = 1.2, y = Inf, label = paste("K-W p =", p.value)), 
            vjust = 1.5, hjust = 0.5, size = 3, inherit.aes = FALSE, color = "steelblue4")
