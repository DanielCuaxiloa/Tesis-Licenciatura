
#########################################################################
# Undirected Gaussian Graphical Model - Quadratic Discriminant Analysis #
# UGGM_QDA                                                              #
#########################################################################

UGGM_QDA <- function(formula, data, rho, prior = NULL) {
  
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

predict.UGGM_QDA <- function(object, newdata){
  
  NewData <- as.matrix(newdata)
  mu <- object$mu
  pi <- object$pi
  Sigma <- object$sigma
  Omega <- object$omega
  
  # Número de clases y observaciones
  num.clases <- length(mu)
  num.obs <- nrow(NewData)
  
  # Probabilidades de pertenecer a cada clase para cada observación
  probs <- matrix(NA, nrow = num.obs, ncol = num.clases)
  
  for (j in 1:num.clases) {
    
    # Función discriminante para la clase j
    delta <- log(pi[[j]]) - 0.5 * log(det(solve(Omega[[j]]))) - 0.5 * rowSums((NewData - rep(mu[[j]], each = num.obs)) %*% Omega[[j]] * (NewData - rep(mu[[j]], each = num.obs)))
    
    # Matriz de probabilidades
    probs[, j] <- delta
    
  }
  
  # Probabilidades normalizadas
  probs <- exp(probs)
  probs <- data.frame(probs = probs/rowSums(probs))
  
  # Nombre de las clases
  nombres.clases <- names(pi)
  
  # Clase predicha para cada observación con los nombres de las clases
  clase.predicha <- data.frame(yhat = apply(probs, 1, function(x) nombres.clases[which.max(x)]))
  
  colnames(probs) <- nombres.clases
  
  # Probabilidad estimada y clase predicha
  return(list(prob = probs, 
              clase = clase.predicha))
  
}

tune.rho <- function(formula, data, rhos, prior = NULL, nfolds = 5) {
  
  etiqueta <- as.character(formula[[2]])
  
  # Resultados de la validación cruzada
  cv.results <- data.frame(rho = rhos, 
                           accuracy = rep(NA, length(rhos)),
                           std.dev = rep(NA, length(rhos)))
  
  # Subconjuntos (folds) de entrenamiento y validación
  set.seed(123)
  folds <- caret::createFolds(data[[etiqueta]], k = nfolds)
  
  for (i in 1:length(rhos)) {
    
    accuracy.val <- numeric(nfolds)
    
    # Realizar validación cruzada
    for (fold_id in folds) {
      
      MC <- NULL
      
      # Conjunto de entrenamiento
      train <- data[-fold_id, ]
      
      # Conjunto de validación
      val <- data[fold_id, ]
      
      # Modelo
      model <- UGGM_QDA(formula, prior = prior, train, rhos[i])
      
      # Predicciones
      predictions <- predict.UGGM_QDA(model, val[, -which(names(data) == etiqueta)])
      
      MC <- table(val[[etiqueta]], predictions$clase$yhat)
      
      # Precisión
      accuracy.val[fold_id] <- sum(diag(MC))/sum(MC)
      
    }
    
    # Precisión promedio
    avg.accuracy <- mean(accuracy.val)
    
    # Desviación estándar de la precisión
    std.dev <- sd(accuracy.val)
    
    # Precisión promedio
    cv.results$accuracy[i] <- avg.accuracy
    cv.results$std.dev[i] <- std.dev
  }
  
  # Valor de rho que maximiza la precisión
  best.rho <- cv.results$rho[which.max(cv.results$accuracy)]
  
  # Resultados
  return(list(cv.results = cv.results, best.rho = best.rho))
}

predict.theoricQDA <- function(object, newdata){
  
  NewData <- as.matrix(newdata)
  mu <- object$mu
  pi <- object$pi
  Sigma <- object$sigma

  num.clases <- length(mu)
  num.obs <- nrow(NewData)
  
  probs <- matrix(NA, nrow = num.obs, ncol = num.clases)
  
  for (j in 1:num.clases) {
    
    delta <- log(pi[[j]]) - 0.5 * log(det(Sigma[[j]])) - 0.5 * rowSums((NewData - rep(mu[[j]], each = num.obs)) %*% solve(Sigma[[j]]) * (NewData - rep(mu[[j]], each = num.obs)))
    
    probs[, j] <- delta
    
  }
  
  probs <- exp(probs)
  probs <- data.frame(probs = probs/rowSums(probs))
  
  nombres.clases <- names(pi)
  
  clase.pred <- data.frame(yhat = apply(probs, 1, function(x) nombres.clases[which.max(x)]))
  
  colnames(probs) <- nombres.clases
  
  return(list(prob = probs, 
              clase = clase.pred))
  
}