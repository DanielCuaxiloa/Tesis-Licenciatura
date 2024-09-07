
###########################
# Simulación 0 | UGGM_QDA #
###########################

library(igraph)
library(qgraph)
library(bootnet)
library(MASS)
library(dplyr)
library(glasso)

library(ggplot2)
library(ggcorrplot)

# RAND --------------------------------------------------------------------

set.seed(12345)

RAND <- function(alpha, p) {
  
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
    diag(adjacency.matrix)[i] <- 1.01*sum(abs(adjacency.matrix[,i]))
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

Modelo <- RAND(alpha = 1.75, p = 10)

plot(x = Modelo$network,
     layout = layout.circle,
     vertex.size = 15,
     vertex.color = "skyblue",
     vertex.shape = "circle",
     edge.color = "gray",
     edge.width = 2,
     main = "Gráfica asociada al modelo UGGM")

aux <- round(Modelo$K,2)
rho <- matrix(0, nrow = nrow(aux), ncol = ncol(aux))

for (i in 1:nrow(aux)) {
  for (j in 1:ncol(aux)) {
    rho[i,j] <- -aux[i,j]/sqrt(aux[i,i]*aux[j,j])
  }
}

Datos_Simulados <- mvrnorm(n = 1000, mu = Modelo$m, Sigma = Modelo$S) %>%
  data.frame()


# Gráficas ----------------------------------------------------------------

cor_matrix <- cor(Datos_Simulados)

ggcorrplot(cor_matrix, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           type = "lower", 
           lab = TRUE, 
           title = "Matriz de correlaciones") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme_bw()

# -------------------------------------------------------------------------

layout(matrix(c(1,2,3,
                4,5,6), 2, 3, byrow = TRUE))

# rho1 = 0.84
# rho2 = 0.86
# rho3 = 0.88
# rho4 = 0.90
# rho5 = 0.92
# rho6 = 0.94

UGGM <- glasso(s = cor_matrix, rho = 0.94)

concentration_matrix <- UGGM$wi

S <- matrix(0, nrow = nrow(concentration_matrix), ncol = ncol(concentration_matrix))

for (i in 1:nrow(S)) {
  for (j in 1:ncol(S)) {
    S[i,j] <- -concentration_matrix[i,j]/sqrt(concentration_matrix[i,i]*concentration_matrix[j,j])
  }
}

S[lower.tri(S)] <- 0

qgraph(S, directed = FALSE, edge.label.cex = 0.6)
title("lambda = 0.94", line = 2)


