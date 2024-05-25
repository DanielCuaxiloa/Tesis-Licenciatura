
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

RAND.1 <- RAND(alpha = 1, p = 30)
RAND.2 <- RAND(alpha = 1.1, p = 30)
RAND.3 <- RAND(alpha = 1.2, p = 30)

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

## n = 50

muestra1.Train1 <- mvrnorm(n = 50, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Train1 <- mvrnorm(n = 50, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Train1 <- mvrnorm(n = 50, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))

Datos.Train1 <- bind_rows(muestra1.Train1, 
                          muestra2.Train1,
                          muestra3.Train1) 

## n = 100

muestra1.Train2 <- mvrnorm(n = 100, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Train2 <- mvrnorm(n = 100, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Train2 <- mvrnorm(n = 100, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))

Datos.Train2 <- bind_rows(muestra1.Train2, 
                          muestra2.Train2,
                          muestra3.Train2)

## n = 1000

muestra1.Train3 <- mvrnorm(n = 1000, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Train3 <- mvrnorm(n = 1000, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Train3 <- mvrnorm(n = 1000, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))

Datos.Train3 <- bind_rows(muestra1.Train3, 
                          muestra2.Train3,
                          muestra3.Train3) 


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

rho.tune1 <- tune.rho(formula = Clase~.,
                      data = Datos.Train1,
                      rhos = seq(0, 1, by = 0.01),
                      nfolds = 5)

rho.tune2 <- tune.rho(formula = Clase~.,
                      data = Datos.Train2,
                      rhos = seq(0, 1, by = 0.01),
                      nfolds = 5)

rho.tune3 <- tune.rho(formula = Clase~.,
                      data = Datos.Train3,
                      rhos = seq(0, 1, by = 0.01),
                      nfolds = 5)

plot(rho.tune1[["cv.results"]]$rho,
     rho.tune1[["cv.results"]]$accuracy)

plot(rho.tune2[["cv.results"]]$rho,
     rho.tune2[["cv.results"]]$accuracy)

plot(rho.tune3[["cv.results"]]$rho,
     rho.tune3[["cv.results"]]$accuracy)

rho.tune1$best.rho
rho.tune2$best.rho
rho.tune3$best.rho


Modelo1.NetQDA <- NetQDA(formula = Clase~.,
                         data = Datos.Train1,
                         rho = rho.tune1$best.rho)

Modelo2.NetQDA <- NetQDA(formula = Clase~.,
                         data = Datos.Train2,
                         rho = rho.tune2$best.rho)

Modelo3.NetQDA <- NetQDA(formula = Clase~.,
                         data = Datos.Train3,
                         rho = rho.tune3$best.rho)

Pred1.NetQDA <- predict.NetQDA(object = Modelo1.NetQDA,
                               newdata = select(Datos.Test, -Clase))

Pred2.NetQDA <- predict.NetQDA(object = Modelo2.NetQDA,
                               newdata = select(Datos.Test, -Clase))

Pred3.NetQDA <- predict.NetQDA(object = Modelo3.NetQDA,
                               newdata = select(Datos.Test, -Clase))

Modelo1.QDA <- qda(formula = Clase~.,
                   data = Datos.Train1)

Modelo2.QDA <- qda(formula = Clase~.,
                   data = Datos.Train2)

Modelo3.QDA <- qda(formula = Clase~.,
                   data = Datos.Train3)

Pred1.QDA <- predict(object = Modelo1.QDA,
                    newdata = select(Datos.Test, -Clase))$class

Pred2.QDA <- predict(object = Modelo2.QDA,
                     newdata = select(Datos.Test, -Clase))$class

Pred3.QDA <- predict(object = Modelo3.QDA,
                     newdata = select(Datos.Test, -Clase))$class

MC1.NetQDA <- table(Datos.Test$Clase, Pred1.NetQDA$clase$yhat)
MC1.QDA <- table(Datos.Test$Clase, Pred1.QDA)

MC2.NetQDA <- table(Datos.Test$Clase, Pred2.NetQDA$clase$yhat)
MC2.QDA <- table(Datos.Test$Clase, Pred2.QDA)

MC3.NetQDA <- table(Datos.Test$Clase, Pred3.NetQDA$clase$yhat)
MC3.QDA <- table(Datos.Test$Clase, Pred3.QDA)


sum(diag(MC1.NetQDA))/sum(MC1.NetQDA)
sum(diag(MC1.QDA))/sum(MC1.QDA)

sum(diag(MC2.NetQDA))/sum(MC2.NetQDA)
sum(diag(MC2.QDA))/sum(MC2.QDA)

sum(diag(MC3.NetQDA))/sum(MC3.NetQDA)
sum(diag(MC3.QDA))/sum(MC3.QDA)


Resultados <- data.frame(
  Modelo = factor(c("NetQDA", "NetQDA", "NetQDA", "QDA", "QDA", "QDA")),
  Rho = c(rho.tune1$best.rho, rho.tune2$best.rho, rho.tune3$best.rho, NA, NA, NA),
  N = c(50, 100, 1000, 50, 100, 1000),
  Precision = c(round(sum(diag(MC1.NetQDA))/sum(MC1.NetQDA),2),
                round(sum(diag(MC2.NetQDA))/sum(MC2.NetQDA),2),
                round(sum(diag(MC3.NetQDA))/sum(MC3.NetQDA),2),
                round(sum(diag(MC1.QDA))/sum(MC1.QDA),2),
                round(sum(diag(MC2.QDA))/sum(MC2.QDA),2),
                round(sum(diag(MC3.QDA))/sum(MC3.QDA),2))
)

ggplot(data = Resultados,
       mapping = aes(x = N, y = Precision, color = Modelo)) +
  geom_line()

# Plot 1 ------------------------------------------------------------------

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
plot(RAND.1$network)
plot(RAND.2$network)
plot(RAND.3$network)


# Plot 2 ------------------------------------------------------------------

network.estimate.Clase1 <- estimateNetwork(data = select(muestra1.Train3,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

network.estimate.Clase2 <- estimateNetwork(data = select(muestra2.Train3,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

network.estimate.Clase3 <- estimateNetwork(data = select(muestra3.Train3,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
plot(network.estimate.Clase1)
plot(network.estimate.Clase2)
plot(network.estimate.Clase3)


