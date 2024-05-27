
###############
# Network QDA #
###############

library(igraph)
library(qgraph)
library(bootnet)
library(MASS)
library(dplyr)


# Estructuras -------------------------------------------------------------

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


# -------------------------------------------------------------------------

Data <- list(Train = list(), 
             Test = list())

Network <- list()

N <- c(50, 100, 500, 1000)
P <- c(20, 25, 30, 35)

for (p in P) {
  
  RAND.1 <- RAND(alpha = 1, p = p)
  RAND.2 <- RAND(alpha = 1.1, p = p)
  RAND.3 <- RAND(alpha = 1.2, p = p)
  
  S1 <- RAND.1$S
  K1 <- solve(S1)
  m1 <- rep(x = 0, times = ncol(S1))
  
  S2 <- RAND.2$S
  K2 <- solve(S2)
  m2 <- rep(x = 0, times = ncol(S2))
  
  S3 <- RAND.3$S
  K3 <- solve(S3)
  m3 <- rep(x = 0, times = ncol(S3))
  
  Network[[as.character(p)]] <- list(RAND.1$network, RAND.2$network, RAND.3$network)
  
  set.seed(123)
  data.train_p <- lapply(N, function(n) {
    
    muestra1 <- mvrnorm(n = n, mu = m1, Sigma = S1) %>% 
      data.frame() %>% 
      mutate(Clase = as.factor("Clase1"))
    
    muestra2 <- mvrnorm(n = n, mu = m2, Sigma = S2) %>% 
      data.frame() %>% 
      mutate(Clase = as.factor("Clase2"))
    
    muestra3 <- mvrnorm(n = n, mu = m3, Sigma = S3) %>% 
      data.frame() %>% 
      mutate(Clase = as.factor("Clase3"))
    
    bind_rows(muestra1, muestra2, muestra3) 
  })
  
  names(data.train_p) <- as.character(N)
  Data$Train[[as.character(p)]] <- data.train_p
  
  set.seed(321)
  
  muestra1_test <- mvrnorm(n = 1000, mu = m1, Sigma = S1) %>% 
    data.frame() %>% 
    mutate(Clase = as.factor("Clase1"))
  
  muestra2_test <- mvrnorm(n = 1000, mu = m2, Sigma = S2) %>% 
    data.frame() %>% 
    mutate(Clase = as.factor("Clase2"))
  
  muestra3_test <- mvrnorm(n = 1000, mu = m3, Sigma = S3) %>% 
    data.frame() %>% 
    mutate(Clase = as.factor("Clase3"))
  
  Data$Test[[as.character(p)]] <- bind_rows(muestra1_test, muestra2_test, muestra3_test)
}

# Funciones NetQDA --------------------------------------------------------

source("NetQDA.R")

# -------------------------------------------------------------------------

rho.tune <- list()

for (p in P) {
  rho.tune_p <- lapply(Data$Train[[as.character(p)]], function(Datos.Train) {
    tune.rho(formula = Clase~.,
             data = Datos.Train,
             rhos = seq(0, 1, by = 0.01),
             nfolds = 10)
  })
  
  names(rho.tune_p) <- as.character(N)
  rho.tune[[as.character(p)]] <- rho.tune_p
}


# -------------------------------------------------------------------------

cv.results <- data.frame()

for (p in P) {
  for (n in N) {
    results <- rho.tune[[as.character(p)]][[as.character(n)]][["cv.results"]]
    results_df <- data.frame(rho = results$rho, 
                             accuracy = results$accuracy, 
                             std.dev = results$std.dev)
    results_df$P <- p
    results_df$N <- n
    cv.results <- bind_rows(cv.results, results_df)
  }
}

best.rho <- cv.results %>% 
  group_by(P, N) %>% 
  summarize(best.rho = rho[which.max(accuracy)], .groups = 'drop')

ggplot(data = cv.results, 
       mapping = aes(x = rho, 
                     y = accuracy)) +
  geom_point(size = 0.1) +
  geom_ribbon(aes(ymin = accuracy - std.dev, 
                  ymax = accuracy + std.dev), 
              alpha = 0.5,
              fill = "grey") +  
  geom_vline(data = best.rho, 
             aes(xintercept = best.rho), 
             linetype = "dashed", 
             color = "red") + 
  facet_grid(N ~ P) +
  theme_bw() +
  labs(title = "Tuning Results",
       x = "Rho",
       y = "Accuracy",
       caption = "Each panel represents a different combination of P and N")


# -------------------------------------------------------------------------

NetQDA.models <- list()
Pred.NetQDA <- list()
MC.NetQDA <- list()
Accuracy.NetQDA <- list()

for (p in P) {
  
  NetQDA.models_p <- list()
  Pred.NetQDA_p <- list()
  MC.NetQDA_p <- list()
  Accuracy.NetQDA_p <- list()
  
  for (n in names(Data$Train[[as.character(p)]])) {
    
    Datos.Train <- Data$Train[[as.character(p)]][[n]]
    best_rho_value <- best.rho %>% filter(P == p & N == as.numeric(n)) %>% pull(best.rho)
    
    model <- NetQDA(formula = Clase~.,
                    data = Datos.Train,
                    rho = best_rho_value)
    
    NetQDA.models_p[[n]] <- model
    
    Datos.Test <- Data$Test[[as.character(p)]]
    predicciones <- predict.NetQDA(object = model, 
                                   newdata = select(Datos.Test, -Clase))
    
    Pred.NetQDA_p[[n]] <- predicciones
    
    confusion_matrix <- table(Datos.Test$Clase, predicciones$clase$yhat)
    MC.NetQDA_p[[n]] <- confusion_matrix
    
    accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
    Accuracy.NetQDA_p[[n]] <- accuracy
  }
  
  #names(NetQDA.models_p) <- names(Data$Train[[as.character(p)]])
  #names(Pred.NetQDA_p) <- names(Data$Train[[as.character(p)]])
  
  NetQDA.models[[as.character(p)]] <- NetQDA.models_p
  Pred.NetQDA[[as.character(p)]] <- Pred.NetQDA_p
  MC.NetQDA[[as.character(p)]] <- MC.NetQDA_p
  Accuracy.NetQDA[[as.character(p)]] <- Accuracy.NetQDA_p
  
}

Resultados.NetQDA <- data.frame(P = rep(c(20, 25, 30, 35), each = 4),
                                N = rep(c(50, 100, 500, 1000), times = 4),
                                Accuracy = unlist(Accuracy.NetQDA),
                                Modelo = factor("NetQDA"))


# -------------------------------------------------------------------------


QDA.models <- list()
Pred.QDA <- list()
MC.QDA <- list()
Accuracy.QDA <- list()

for (p in P) {
  
  QDA.models_p <- list()
  Pred.QDA_p <- list()
  MC.QDA_p <- list()
  Accuracy.QDA_p <- list()
  
  for (n in names(Data$Train[[as.character(p)]])) {
    
    Datos.Train <- Data$Train[[as.character(p)]][[n]]

    model <- qda(formula = Clase~.,
                 data = Datos.Train)
    
    QDA.models_p[[n]] <- model
    
    Datos.Test <- Data$Test[[as.character(p)]]
    predicciones <- predict(object = model, 
                            newdata = select(Datos.Test, -Clase))$class
    
    Pred.QDA_p[[n]] <- predicciones
    
    confusion_matrix <- table(Datos.Test$Clase, predicciones)
    MC.QDA_p[[n]] <- confusion_matrix
    
    accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
    Accuracy.QDA_p[[n]] <- accuracy
  }
  
  #names(NetQDA.models_p) <- names(Data$Train[[as.character(p)]])
  #names(Pred.NetQDA_p) <- names(Data$Train[[as.character(p)]])
  
  QDA.models[[as.character(p)]] <- QDA.models_p
  Pred.QDA[[as.character(p)]] <- Pred.QDA_p
  MC.QDA[[as.character(p)]] <- MC.QDA_p
  Accuracy.QDA[[as.character(p)]] <- Accuracy.QDA_p
  
}

Resultados.QDA <- data.frame(P = rep(c(20, 25, 30, 35), each = 4),
                             N = rep(c(50, 100, 500, 1000), times = 4),
                             Accuracy = unlist(Accuracy.QDA),
                             Modelo = factor("QDA"))

# -------------------------------------------------------------------------

Resultados <- bind_rows(Resultados.NetQDA,
                        Resultados.QDA) %>% 
  mutate(Accuracy = round(Accuracy,4))

rownames(Resultados) <- NULL

ggplot(data = Resultados,
       mapping = aes(x = as.factor(N), 
                     y = Accuracy, 
                     color = Modelo,
                     shape = Modelo)) +
  geom_point() +
  geom_line(aes(group = Modelo)) + 
  facet_grid(~P) + 
  theme_bw() + 
  labs(x = "n", 
       y = "Precisión", 
       color = "Modelo", 
       shape = "Modelo")



