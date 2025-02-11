
###########################
# Simulación 1 | UGGM_QDA #
###########################

library(igraph)
library(qgraph)
library(bootnet)
library(MASS)
library(dplyr)
library(glasso)

library(ggplot2)
library(ggcorrplot)
library(ppcor)


# RAND --------------------------------------------------------------------

set.seed(12345)

RAND <- function(alpha, p, R) {
  
  network <- sample_pa(n = p, power = alpha, directed = FALSE)
  
  adjacency.matrix <- as.matrix(as_adjacency_matrix(network))
  
  aux <- function() {
    if (runif(1) < 0.5) {
      return(runif(1, -1, -0.5))
    } else {
      return(runif(1, 0.5, 1))
    }
  }
  
  for(i in 1:nrow(adjacency.matrix)) {
    for(j in 1:i) {
      if(adjacency.matrix[i, j] != 0){
        adjacency.matrix[i, j] <- adjacency.matrix[j, i] <- aux()
      } else{
        adjacency.matrix[i, j] <- adjacency.matrix[i, j] 
      }
    }
  }
  
  for (i in 1:p) {
    diag(adjacency.matrix)[i] <- R*sum(abs(adjacency.matrix[,i]))
  }
  
  Aux <- solve(adjacency.matrix)
  
  S <- matrix(0, nrow = nrow(Aux), ncol = ncol(Aux))
  
  for (i in 1:nrow(S)) {
    for (j in 1:ncol(S)) {
      S[i,j] <- Aux[i,j]/sqrt(Aux[i,i]*Aux[j,j])
    }
  }
  
  K <- solve(S)
  m <- rep(0, ncol(S))
  
  return(list(S = S, K = K, m = m, network = network))
  
}

Modelo <- RAND(alpha = 1.75, p = 10, R = 1.01)

node_labels <- paste0("v[", 1:vcount(Modelo$network), "]")
node_labels <- parse(text = node_labels)

plot(x = Modelo$network,
     layout = layout_in_circle(Modelo$network), 
     vertex.size = 25,                   
     vertex.color = "skyblue",           
     vertex.frame.color = "black",       
     vertex.label.color = "black",       
     vertex.label = node_labels,
     vertex.label.cex = 1.5,             
     vertex.label.font = 2,              
     edge.color = "gray50",              
     edge.width = 2                      
)

Datos_Simulados <- mvrnorm(n = 1000, mu = Modelo$m, Sigma = Modelo$S) %>%
  data.frame()


# Gráficas ----------------------------------------------------------------

cor_matrix <- pcor(Datos_Simulados)

ggcorrplot(cor_matrix$estimate, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           type = "lower", 
           lab = TRUE, 
           lab_size = 5,
           lab_col = "black") +
  labs(title = "Matriz de Correlaciones Parciales de Pearson",
       subtitle = "Datos Simulados",
       x = "",
       y = "") + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

qgraph(cor_matrix$estimate,
       layout="circle",
       layout = "spring", 
       labels = colnames(Datos_Simulados), 
       theme = "classic",
       edge.color = ifelse(cor_matrix$estimate > 0, "#E46726", "#6D9EC1"),
       title = "Modelo Gráfico basado en Correlaciones Parciales",
       title.cex = 1.5)


# glasso ------------------------------------------------------------------

layout(matrix(c(1, 2,
                3, 4,
                5, 6), 3, 2, byrow = TRUE))

# rho1 = 0.84
# rho2 = 0.86
# rho3 = 0.88
# rho4 = 0.90
# rho5 = 0.92
# rho6 = 0.94

UGGM <- glasso(s = cov(Datos_Simulados), rho = 0.94)

concentration_matrix <- UGGM$wi

S <- matrix(0, 
            nrow = nrow(concentration_matrix), 
            ncol = ncol(concentration_matrix))

for (i in 1:nrow(S)) {
  for (j in 1:ncol(S)) {
    S[i,j] <- -concentration_matrix[i,j]/sqrt(concentration_matrix[i,i]*concentration_matrix[j,j])
  }
}

S[lower.tri(S)] <- 0

qgraph(S,
       directed = FALSE,
       layout="circle",
       layout = "spring", 
       labels = colnames(Datos_Simulados), 
       theme = "classic",
       edge.color = ifelse(cor_matrix$estimate > 0, "#E46726", "#6D9EC1"),
       title = bquote(lambda == 0.94),
       title.cex = 1.5,
       label.cex = 1.5)



