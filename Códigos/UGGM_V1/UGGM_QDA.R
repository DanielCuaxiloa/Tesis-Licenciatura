
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
  
  # Verificación de que el vector de rho tiene el tamaño adecuado
  if (length(rho) != length(niveles.y)) {
    stop("El número de hiperparámetros rho debe coincidir con el número de clases en la variable y")
  }
  
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
  
  # Matriz de precisión por cada clase utilizando el valor de rho correspondiente
  Omega <- mapply(function(sigma, rho_nivel) {
    glasso::glasso(sigma, rho = rho_nivel)$wi
  }, Sigma, rho, SIMPLIFY = FALSE)
  
  # Asignar nombres de las clases
  names(mu) <- niveles.y
  names(pi) <- niveles.y
  names(Sigma) <- niveles.y
  names(Omega) <- niveles.y
  
  # Estimaciones
  return(list(mu = mu, pi = pi, sigma = Sigma, omega = Omega))
}

predict.UGGM_QDA <- function(object, newdata) {
  
  # Convertir newdata a matriz
  NewData <- as.matrix(newdata)
  
  # Extraer parámetros del modelo entrenado
  mu <- object$mu
  pi <- object$pi
  Omega <- object$omega
  
  # Número de clases y observaciones
  num.clases <- length(mu)
  num.obs <- nrow(NewData)
  
  # Matriz para almacenar las probabilidades de pertenecer a cada clase
  probs <- matrix(NA, nrow = num.obs, ncol = num.clases)
  
  # Calcular la función discriminante para cada clase
  for (j in 1:num.clases) {
    
    # Función discriminante para la clase j
    delta <- log(pi[[j]]) - 0.5 * log(det(Omega[[j]])) - 
      0.5 * rowSums((NewData - matrix(mu[[j]], nrow = num.obs, ncol = length(mu[[j]]), byrow = TRUE)) %*% Omega[[j]] * 
                      (NewData - matrix(mu[[j]], nrow = num.obs, ncol = length(mu[[j]]), byrow = TRUE)))
    
    # Almacenar las probabilidades discriminantes
    probs[, j] <- delta
  }
  
  # Convertir las probabilidades discriminantes a probabilidades reales
  probs <- exp(probs)
  probs <- probs / rowSums(probs) # Normalización
  
  # Convertir a data.frame para mejor manejo
  probs <- data.frame(probs)
  
  # Nombres de las clases (etiquetas)
  nombres.clases <- names(pi)
  colnames(probs) <- nombres.clases
  
  # Predecir la clase con la mayor probabilidad
  clase.predicha <- apply(probs, 1, function(x) nombres.clases[which.max(x)])
  
  # Convertir a data.frame
  clase.predicha <- data.frame(yhat = clase.predicha)
  
  # Retornar probabilidades y predicciones
  return(list(prob = probs, clase = clase.predicha))
}

tune.rho <- function(formula, data, rhos.list, prior = NULL, nfolds = 5) {
  
  etiqueta <- as.character(formula[[2]])
  niveles.y <- unique(data[[etiqueta]])
  
  # Verificar que haya tantos conjuntos de valores rho como clases
  if (length(rhos.list) != length(niveles.y)) {
    stop("El número de conjuntos de valores de rho debe coincidir con el número de clases.")
  }
  
  # Generar todas las combinaciones posibles de rhos
  combinaciones.rho <- expand.grid(rhos.list)
  
  # Resultados de la validación cruzada
  cv.results <- data.frame(matrix(NA, nrow = nrow(combinaciones.rho), ncol = length(niveles.y) + 2))
  colnames(cv.results) <- c(paste0("rho_clase_", seq_along(niveles.y)), "accuracy", "std.dev")
  
  # Crear subconjuntos (folds) de entrenamiento y validación
  set.seed(123)
  folds <- caret::createFolds(data[[etiqueta]], k = nfolds)
  
  # Validar para cada combinación de valores de rho
  for (i in 1:nrow(combinaciones.rho)) {
    
    # Almacenar precisión en cada fold
    accuracy.val <- numeric(nfolds)
    
    # Obtener la combinación actual de rhos
    rho.combinacion <- as.numeric(combinaciones.rho[i, ])
    
    # Realizar validación cruzada
    for (fold_index in seq_along(folds)) {
      
      fold_id <- folds[[fold_index]]
      
      # Conjunto de entrenamiento
      train <- data[-fold_id, ]
      
      # Conjunto de validación
      val <- data[fold_id, ]
      
      # Entrenar modelo con la combinación actual de rhos
      model <- UGGM_QDA(formula = formula, data = train, rho = rho.combinacion, prior = prior)
      
      # Predicciones sobre el conjunto de validación
      predictions <- predict.UGGM_QDA(model, newdata = val[, -which(names(val) == etiqueta)])
      
      # Matriz de confusión
      MC <- table(val[[etiqueta]], predictions$clase$yhat)
      
      # Calcular la precisión
      accuracy.val[fold_index] <- sum(diag(MC)) / sum(MC)
    }
    
    # Precisión promedio y desviación estándar
    avg.accuracy <- mean(accuracy.val)
    std.dev <- sd(accuracy.val)
    
    # Guardar los resultados en el dataframe
    cv.results[i, 1:length(niveles.y)] <- rho.combinacion
    cv.results$accuracy[i] <- avg.accuracy
    cv.results$std.dev[i] <- std.dev
  }
  
  # Seleccionar la combinación de rhos que maximiza la precisión
  best.rho <- combinaciones.rho[which.max(cv.results$accuracy), ]
  
  # Devolver resultados de la validación cruzada y la mejor combinación de rhos
  return(list(cv.results = cv.results, best.rho = best.rho))
}




