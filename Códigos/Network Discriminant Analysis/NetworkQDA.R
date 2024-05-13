
###############
# Network QDA #
###############

library(igraph)
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
  
  return(S)
  
}


# Simulaciones ------------------------------------------------------------

S1 <- RAND(alpha = 0.5, p = 30)
K1 <- solve(S1)
m1 <- rep(x = 0, times = ncol(S1))

S2 <- RAND(alpha = 1, p = 30)
K2 <- solve(S2)
m2 <- rep(x = 0, times = ncol(S2))

S3 <- RAND(alpha = 1.5, p = 30)
K3 <- solve(S3)
m3 <- rep(x = 0, times = ncol(S3))


# Train -------------------------------------------------------------------

set.seed(123)

muestra1.Train <- mvrnorm(n = 500, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2.Train <- mvrnorm(n = 500, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3.Train <- mvrnorm(n = 500, mu = m3, Sigma = S3) %>% 
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


# Función NetQDA ----------------------------------------------------------

NetQDA <- function(formula, data, rho, prior = NULL) {
  
  etiqueta <- as.character(formula[[2]])
  
  if (!inherits(formula, "formula")) {
    stop("El argumento 'formula' debe ser tipo formula")
  }
  
  if (!is.data.frame(data)) {
    stop("El argumento 'datos' debe ser un data.frame.")
  }
  
  # Si el lado derecho de la fórmula es '.', usar todas las otras columnas como predictores
  if (length(formula[[3]]) == 1 && as.character(formula[[3]]) == ".") {
    predictores <- setdiff(names(data), etiqueta)
    formula <- reformulate(predictores, response = etiqueta)
  }
  
  # Vector de etiquetas
  y <- as.factor(data[[etiqueta]])
  
  # Matriz de diseño sin intercepto
  X <- model.matrix(formula, data)[,-1]
  
  # Niveles únicos de la variable y
  niveles.y <- unique(y)
  
  # Medias de cada variable en X agrupadas por las clases de y
  mu <- lapply(niveles.y, function(nivel) {
    colMeans(X[y == nivel, , drop = FALSE])
  })
  
  # Probabilidades a priori por cada clase
  if (is.null(prior)) {
    pi <- as.list(table(y)/length(y))
  } else {
    if (length(prior) != length(niveles.y)) {
      stop("El número de probabilidades a priori no coincide con el número de clases en la variable de y")
    }
    pi <- as.list(prior/sum(prior))
  }
  
  # Matriz de covarianza por cada clase
  Sigma <- lapply(niveles.y, function(nivel) {
    as.matrix(var(X[y == nivel, , drop = FALSE]))
  })
  
  # Matriz de precisión por cada clase
  Omega <- lapply(Sigma, function(sigma) {
    glasso::glasso(sigma, rho = rho)$wi
  })
  
  # Asignar nombres de las clases
  names(mu) <- niveles.y
  
  names(pi) <- niveles.y
  
  names(Sigma) <- niveles.y
  
  names(Omega) <- niveles.y
  
  # Estimaciones
  return(list(mu = mu, pi = pi, sigma = Sigma, omega = Omega))
}

predict.NetQDA <- function(object, newdata){
  
  NewData <- as.matrix(newdata)
  mu <- object$mu
  pi <- object$pi
  Sigma <- object$sigma
  Omega <- object$omega
  
  # Número de clases y observaciones
  num.clases <- length(mu)
  num.obs <- nrow(NewData)
  
  # Probabilidades de pertenecer a cada clase para cada observación
  probabilidades <- matrix(NA, nrow = num.obs, ncol = num.clases)
  
  for (j in 1:num.clases) {
    
    # Función discriminante para la clase j
    delta <- log(pi[[j]]) - 0.5 * log(det(solve(Omega[[j]]))) - 0.5 * rowSums((NewData - rep(mu[[j]], each = num.obs)) %*% Omega[[j]] * (NewData - rep(mu[[j]], each = num.obs)))
    
    # Matriz de probabilidades
    probabilidades[, j] <- delta
    
    }
  
  # Probabilidades normalizadas
  probabilidades <- exp(probabilidades)
  probabilidades <- data.frame(probabilidades = probabilidades/rowSums(probabilidades))
  
  # Nombre de las clases
  nombres.clases <- names(pi)
  
  # Clase predicha para cada observación con los nombres de las clases
  clase.predicha <- data.frame(yhat = apply(probabilidades, 1, function(x) nombres.clases[which.max(x)]))
  
  colnames(probabilidades) <- nombres.clases
  
  # Probabilidad estimada y clase predicha
  return(list(probabilidad = probabilidades, 
              clase = clase.predicha))
  
}

tune.rho <- function(formula, data, rhos, prior = NULL, nfolds = 5) {
  
  etiqueta <- as.character(formula[[2]])
  
  # Resultados de la validación cruzada
  cv.results <- data.frame(rho = rhos, accuracy = rep(NA, length(rhos)))
  
  # Subconjuntos (folds) de entrenamiento y validación
  set.seed(123)
  folds <- caret::createFolds(data[[etiqueta]], k = nfolds)
  
  for (i in 1:length(rhos)) {
    
    avg.accuracy <- 0
    
    # Realizar validación cruzada
    for (fold in folds) {
      
      MC <- NULL
      
      # Conjunto de entrenamiento
      train <- data[-fold, ]
      
      # Conjunto de prueba
      val <- data[fold, ]
      
      # Modelo
      model <- NetQDA(formula, prior = prior, train, rhos[i])
      
      # Predicciones
      predictions <- predict.NetQDA(model, val[, -which(names(data) == etiqueta)])
      
      MC <- table(val[[etiqueta]], predictions$clase$yhat)
      
      # Precisión
      accuracy <- sum(diag(MC))/sum(MC)
      
      # Precisión acumulada
      avg.accuracy <- avg.accuracy + accuracy
    }
    
    # Precisión promedio
    cv.results$accuracy[i] <- avg.accuracy/nfolds
  }
  
  # Valor de rho que maximiza la precisión
  best.rho <- cv.results$rho[which.max(cv.results$accuracy)]
  
  # Resultados
  return(list(cv.results = cv.results, best.rho = best.rho))
}


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

network.estimate.Clase1 <- estimateNetwork(data = select(muestra1.Train,-Clase),
                                           default = "EBICglasso",
                                           tuning = 0)

plot(network.estimate.Clase1)
