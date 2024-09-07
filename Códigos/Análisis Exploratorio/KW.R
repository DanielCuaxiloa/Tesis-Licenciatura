# Cargar datos
data <- Datos

# Modificar los nombres de las clases para que sean más claros y ordenarlos
data$Clase <- factor(recode(data$Clase,
                            "GTEX_Brain" = "GTEX Brain",
                            "TCGA_LGG" = "TCGA LGG",
                            "TCGA_GM" = "TCGA GM"),
                     levels = c("GTEX Brain", "TCGA LGG", "TCGA GM"))

# Lista de variables a analizar
variables <- names(data)[1:(ncol(data) - 1)]
group_column <- "Clase"

# Inicializar una lista para almacenar los resultados de pairwise de todas las variables
pairwise_results <- list()
kruskal_results <- data.frame()

# Realizar las comparaciones pairwise para cada variable
for (var in variables) {
  # Convertir la fórmula en un formato apropiado
  formula <- as.formula(paste(var, "~", group_column))
  
  # Prueba de Kruskal-Wallis
  kruskal_test <- kruskal_test(data, formula = formula)
  
  # Añadir los resultados a un dataframe para luego usarlos en la gráfica
  kruskal_results <- rbind(kruskal_results, data.frame(variable = var, p.value = kruskal_test$p))
  
  # Comparaciones pairwise usando Wilcoxon (post-hoc)
  pairwise_test <- pairwise_wilcox_test(data, formula = formula, p.adjust.method = "bonferroni")
  
  # Añadir el nombre de la variable para facilitar el facetado
  pairwise_test$variable <- var
  
  # Añadir la columna y.position para posicionar los p-valores
  max_y <- max(data[[var]], na.rm = TRUE)
  step_increase <- 0.1 * max_y
  pairwise_test <- pairwise_test %>%
    mutate(y.position = seq(max_y + step_increase, by = step_increase, length.out = nrow(pairwise_test)),
           group1 = recode(group1, "GTEX_Brain" = "GTEX Brain", "TCGA_LGG" = "TCGA LGG", "TCGA_GM" = "TCGA GM"),
           group2 = recode(group2, "GTEX_Brain" = "GTEX Brain", "TCGA_LGG" = "TCGA LGG", "TCGA_GM" = "TCGA GM"))
  
  # Guardar los resultados en la lista
  pairwise_results[[var]] <- pairwise_test
}

# Combinar todos los resultados en un solo data frame
pairwise_results_combined <- do.call(rbind, pairwise_results)

# Convertir el data frame original a largo para que funcione con facet_wrap
data_long <- data %>%
  pivot_longer(cols = variables, names_to = "variable", values_to = "value")

# Crear el gráfico combinado con facet_wrap, colorear por clase y mostrar el p-valor de Kruskal-Wallis
ggplot(data_long, aes(x = Clase, y = value)) +
  geom_boxplot(aes(fill = Clase)) +
  stat_pvalue_manual(pairwise_results_combined, label = "p.adj", hide.ns = FALSE, size = 3) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Clase", y = "log2(norm count + 1)", title = "Comparaciones Pairwise con Kruskal-Wallis y Wilcoxon") +
  theme_bw() +
  theme(legend.position = "none") +  # Eliminar la leyenda
  geom_text(data = kruskal_results, aes(x = 1, y = Inf, label = paste("K-W p =", p.value)),
            vjust = 1.5, hjust = 0.5, size = 3, inherit.aes = FALSE, color = "steelblue4")





