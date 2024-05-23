
###############
# Network QDA #
###############

library(igraph)
library(qgraph)
library(bootnet)
library(MASS)
library(dplyr)


# Estructuras -------------------------------------------------------------

AR.1 <- function(rho, p){
  
  S <- matrix(data = NA, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      S[i,j] <- rho^(abs(i-j))
    }
  }
  
  return(S)
  
}

MA.1 <- function(rho, p){
  S <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(abs(i-j)<=1){
        S[i,j] <- rho^(abs(i-j))
      } else{
        S[i,j] <- 0
      }
    }
  }
  return(S)
}

RAND <- function(alpha, p) {
  
  set.seed(123)
  
  # Estructura de gráfico
  network <- sample_pa(n = p, power = alpha, m = 1, directed = FALSE)
  
  # Matriz de adyacencia
  adjacency.matrix <- as.matrix(as_adjacency_matrix(network))
  
  # Definir la función para generar números aleatorios dentro del rango especificado
  aux <- function() {
    if (runif(1) < 0.5) {
      return(runif(1, -1, -0.5))
    } else {
      return(runif(1, 0.5, 1))
    }
  }
  
  # Generar números aleatorios en la matriz de adyacencia
  for(i in 1:nrow(adjacency.matrix)) {
    for(j in 1:i) {
      if(adjacency.matrix[i, j] != 0){
        adjacency.matrix[i, j] <- adjacency.matrix[j, i] <- aux()
      } else{
        adjacency.matrix[i, j] <- adjacency.matrix[i, j] 
      }
    }
  }
  
  # Matriz de adyacencia definida positiva
  for (i in 1:p) {
    diag(adjacency.matrix)[i] <- 1.01*sum(abs(adjacency.matrix[,i]))
  }
  
  # Matriz auxiliar
  Aux <- solve(adjacency.matrix)
  
  S <- matrix(0, nrow = nrow(Aux), ncol = ncol(Aux))
  
  for (i in 1:nrow(S)) {
    for (j in 1:ncol(S)) {
      S[i,j] <- Aux[i,j]/sqrt(Aux[i,i]*Aux[j,j])
    }
  }
  
  return(list(S = S,
              network = network))
  
}


# Simulaciones ------------------------------------------------------------

RAND.1 <- RAND(alpha = 1, p = 50)
RAND.2 <- RAND(alpha = 1.1, p = 50)
RAND.3 <- RAND(alpha = 1.2, p = 50)

S1 <- RAND.1$S
K1 <- solve(S1)
m1 <- rep(x = 0, times = ncol(S1))

S2 <- RAND.2$S
K2 <- solve(S2)
m2 <- rep(x = 0, times = ncol(S2))

S3 <- RAND.3$S
K3 <- solve(S3)
m3 <- rep(x = 0, times = ncol(S3))


# Train -------------------------------------------------------------------

set.seed(123)

muestra1.Train <- mvrnorm(n = 100, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Train <- mvrnorm(n = 100, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Train <- mvrnorm(n = 100, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))

Datos.Train <- bind_rows(muestra1.Train, 
                         muestra2.Train,
                         muestra3.Train) 


# Test --------------------------------------------------------------------

set.seed(321)

muestra1.Test <- mvrnorm(n = 1000, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Test <- mvrnorm(n = 1000, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Test <- mvrnorm(n = 1000, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))


Datos.Test <- bind_rows(muestra1.Test, 
                        muestra2.Test,
                        muestra3.Test) 


# Funciones NetQDA --------------------------------------------------------

source("NetQDA.R")


# Comparación MASS v.s NetQDA ----------------------------------------------

rho.tune <- tune.rho(formula = Clase~.,
                     data = Datos.Train,
                     rhos = seq(0, 1, by = 0.01),
                     nfolds = 5)

plot(rho.tune[["cv.results"]]$rho,rho.tune[["cv.results"]]$accuracy)
rho.tune$best.rho

Modelo.NetQDA <- NetQDA(formula = Clase~.,
                        data = Datos.Train,
                        rho = rho.tune$best.rho)

Pred.NetQDA <- predict.NetQDA(object = Modelo.NetQDA,
                              newdata = select(Datos.Test, -Clase))

Modelo.QDA <- qda(formula = Clase~.,
                  data = Datos.Train)

Pred.QDA <- predict(object = Modelo.QDA,
                    newdata = select(Datos.Test, -Clase))$class

MC.NetQDA <- table(Datos.Test$Clase, Pred.NetQDA$clase$yhat)
MC.QDA <- table(Datos.Test$Clase, Pred.QDA)

sum(diag(MC.NetQDA))/sum(MC.NetQDA)
sum(diag(MC.QDA))/sum(MC.QDA)


# Plot 1 ------------------------------------------------------------------

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
plot(RAND.1$network)
plot(RAND.2$network)
plot(RAND.3$network)


# Plot 2 ------------------------------------------------------------------

network.estimate.Clase1 <- estimateNetwork(data = select(muestra1.Train,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

network.estimate.Clase2 <- estimateNetwork(data = select(muestra2.Train,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

network.estimate.Clase3 <- estimateNetwork(data = select(muestra3.Train,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
plot(network.estimate.Clase1)
plot(network.estimate.Clase2)
plot(network.estimate.Clase3)
