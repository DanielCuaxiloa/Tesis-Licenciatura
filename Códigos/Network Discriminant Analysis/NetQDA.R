
###########################################
# Network Quadratic Discriminant Analysis #
# Network QDA                             #
###########################################

NetQDA <- function(formula, datos, rho) {
  
  # Verificar que el argumento 'formula' sea realmente una fórmula
  if (!inherits(formula, "formula")) {
    stop("El argumento 'formula' debe ser del tipo formula.")
  }
  
  # Verificar que 'datos' sea un data.frame
  if (!is.data.frame(datos)) {
    stop("El argumento 'datos' debe ser un data.frame.")
  }
  
  # Extraer la variable de respuesta y los predictores
  respuesta <- as.character(formula[[2]])
  
  # Si el lado derecho de la fórmula es '.', usar todas las otras columnas como predictores
  if (length(formula[[3]]) == 1 && as.character(formula[[3]]) == ".") {
    predictores <- setdiff(names(datos), respuesta)
    formula <- reformulate(predictores, response = respuesta)
  }
  
  # Extraer el vector de respuesta
  y <- as.factor(datos[[respuesta]])
  
  # Construir la matriz de diseño X usando model.matrix para asegurar correcto manejo de factores
  X <- model.matrix(formula, datos)[,-1]  # Excluye la columna intercepto si no se necesita
  
  # Obtener los niveles únicos de la variable de respuesta
  niveles_y <- unique(y)
  
  # Calcular las medias de cada variable en X agrupadas por las clases de Y
  mu <- lapply(niveles_y, function(nivel) {
    colMeans(X[y == nivel, , drop = FALSE])
  })
  
  # Asignar nombres a las listas basados en los niveles de la variable de respuesta
  names(mu) <- niveles_y
  
  # Calcular las probabilidades a priori por cada clase
  pi <- as.list(table(y)/length(y))
  
  # Asignar nombres a las listas basados en los niveles de la variable de respuesta
  names(pi) <- niveles_y
  
  # Calcular las matrices de covarianza por cada clase
  Sigma <- lapply(niveles_y, function(nivel) {
    as.matrix(var(X[y == nivel, , drop = FALSE]))
  })
  
  names(Sigma) <- niveles_y
  
  # Aplicar glasso() a cada matriz de covarianza
  Omega <- lapply(Sigma, function(sigma) {
    glasso::glasso(sigma, rho = rho)$wi
  })
  
  # Asignar nombres a las listas basados en los niveles de la variable de respuesta
  names(Omega) <- niveles_y
  
  # Devolver los resultados como una lista de medias con nombres
  return(list(y = y, X = X, 
              mu = mu, 
              pi = pi, 
              sigma = Sigma,
              omega = Omega))
}

Predict.NetQDA <- function(object, NewData){
  
  NewData <- as.matrix(NewData)
  mu <- object$mu
  pi <- object$pi
  Sigma <- object$sigma
  Omega <- object$omega
  
  # Calcular el número de clases y observaciones
  num_clases <- length(mu)
  num_obs <- nrow(NewData)
  
  # Calcular las probabilidades de pertenecer a cada clase para cada observación
  probabilidades <- matrix(NA, nrow = num_obs, ncol = num_clases)
  for (j in 1:num_clases) {
    # Calcular la función discriminante para la clase j
    delta <- log(pi[[j]]) - 0.5 * log(det(solve(Omega[[j]]))) - 0.5 * rowSums((NewData - rep(mu[[j]], each = num_obs)) %*% Omega[[j]] * (NewData - rep(mu[[j]], each = num_obs)))
    
    # Almacenar la función discriminante en la matriz de probabilidades
    probabilidades[, j] <- delta
  }
  
  # Calcular las probabilidades normalizadas
  probabilidades <- exp(probabilidades)
  probabilidades <- probabilidades / rowSums(probabilidades)
  
  # Obtener el nombre de las clases
  nombres_clases <- names(pi)
  
  # Obtener la clase predicha para cada observación con los nombres de las clases
  clase_predicha <- apply(probabilidades, 1, function(x) nombres_clases[which.max(x)]) %>% 
    as.data.frame()
  
  colnames(probabilidades) <- nombres_clases
  
  # Devolver la probabilidad estimada y la clase predicha
  return(list(probabilidad = probabilidades, 
              clase_predicha = clase_predicha))
  
}
