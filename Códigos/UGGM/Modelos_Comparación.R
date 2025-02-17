
##########################
# Comparación de Modelos #
##########################

library(ranger)
library(e1071)
library(MASS)
library(glmnet)

library(dplyr)
library(ggplot2)

load("Datos.RData")
source("Modelo_UGGM-QDA.R")

Aux <- function(model.func, predict.func, data.train, data.test, best.rho = NULL) {
  
  models <- list()
  pred <- list()
  MC <- list()
  accuracies <- list()
  
  for (n in names(data.train)) {
    
    Datos.Train <- data.train[[n]]
    
    if (!is.null(best.rho)) {
      best.rho.value <- best.rho %>% filter(P == p & N == as.numeric(n)) %>% pull(best.rho)
      model <- model.func(formula = Clase~., 
                          data = Datos.Train, 
                          rho = best.rho.value)
    } else {
      model <- model.func(formula = Clase~., 
                          data = Datos.Train)
    }
    
    models[[n]] <- model
    
    predicciones <- predict.func(model, select(data.test, -Clase))
    
    pred[[n]] <- predicciones
    
    if (!is.null(best.rho)) {
      mc <- table(data.test$Clase, predicciones$clase$yhat)
    } else {
      mc <- table(data.test$Clase, predicciones)
    }
    
    MC[[n]] <- mc
    accuracies[[n]] <- sum(diag(mc)) / sum(mc)
  }
  
  list(models = models, pred = pred, MC = MC, accuracies = accuracies)
}

Models <- list(UGGM_QDA = list(Model = list(), 
                               Pred = list(), 
                               MC = list(), 
                               Accuracy = list()),
               QDA = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()),
               RF = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()),
               SVM = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()),
               LDA = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()),
               MNL = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()))

rho.tune <- list()

for (p in P) {
  
  rho.tune[[as.character(p)]] <- lapply(gen.data$Data$Train[[as.character(p)]], 
                                        function(Datos.Train) {
                                          tune.rho(formula = Clase~.,
                                                   data = Datos.Train,
                                                   rhos = seq(0, 0.1, by = 0.001),
                                                   nfolds = 10)
                                        })
}

cv.results <- do.call(bind_rows, lapply(names(rho.tune), function(p) {
  do.call(bind_rows, lapply(names(rho.tune[[p]]), function(n) {
    results <- rho.tune[[p]][[n]][["cv.results"]]
    data.frame(P = p,
               N = n,
               rho = results$rho,
               accuracy = results$accuracy,
               std.dev = results$std.dev)
  }))
}))

cv.results$N <- factor(cv.results$N, levels=c("50", "100", "500", "1000"))
cv.results$P <- factor(cv.results$P, levels=c("20", "25", "30", "35"))

best.rho <- cv.results %>% 
  group_by(P, N) %>% 
  summarize(best.rho = rho[which.max(accuracy)], .groups = 'drop')

# UGGM-QDA
for (p in P) {
  result <- Aux(
    model.func = UGGM_QDA,
    predict.func = function(model, newdata) predict.UGGM_QDA(object = model, newdata = newdata),
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]],
    best.rho = best.rho
  )
  
  Models$UGGM_QDA$Model[[as.character(p)]] <- result$models
  Models$UGGM_QDA$Pred[[as.character(p)]] <- result$pred
  Models$UGGM_QDA$MC[[as.character(p)]] <- result$MC
  Models$UGGM_QDA$Accuracy[[as.character(p)]] <- result$accuracies
}

# QDA
for (p in P) {
  result <- Aux(
    model.func = qda,
    predict.func = function(model, newdata) predict(object = model, newdata = newdata)$class,
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]]
  )
  
  Models$QDA$Model[[as.character(p)]] <- result$models
  Models$QDA$Pred[[as.character(p)]] <- result$pred
  Models$QDA$MC[[as.character(p)]] <- result$MC
  Models$QDA$Accuracy[[as.character(p)]] <- result$accuracies
}

# RF
for (p in P) {
  result <- Aux(
    model.func = function(formula, data) ranger(formula, data, probability = FALSE),
    predict.func = function(model, newdata) predict(model, data = newdata)$predictions,
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]]
  )
  
  Models$RF$Model[[as.character(p)]] <- result$models
  Models$RF$Pred[[as.character(p)]] <- result$pred
  Models$RF$MC[[as.character(p)]] <- result$MC
  Models$RF$Accuracy[[as.character(p)]] <- result$accuracies
}

# SVM
for (p in P) {
  result <- Aux(
    model.func = svm,
    predict.func = function(model, newdata) predict(model, newdata),
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]]
  )
  
  Models$SVM$Model[[as.character(p)]] <- result$models
  Models$SVM$Pred[[as.character(p)]] <- result$pred
  Models$SVM$MC[[as.character(p)]] <- result$MC
  Models$SVM$Accuracy[[as.character(p)]] <- result$accuracies
}

# LDA
for (p in P) {
  result <- Aux(
    model.func = lda,
    predict.func = function(model, newdata) predict(model, newdata)$class,
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]]
  )
  
  Models$LDA$Model[[as.character(p)]] <- result$models
  Models$LDA$Pred[[as.character(p)]] <- result$pred
  Models$LDA$MC[[as.character(p)]] <- result$MC
  Models$LDA$Accuracy[[as.character(p)]] <- result$accuracies
}

# MLR (GLMNET con lambda = 0)
for (p in P) {
  result <- Aux(
    model.func = function(formula, data) {
      if (!"Clase" %in% colnames(data)) {
        stop("Error: La variable 'Clase' no está en los datos de entrenamiento")
      }
      data$Clase <- as.factor(data$Clase)
      # Guardar los niveles de Clase
      clase_levels <- levels(data$Clase)
      # Obtener términos del modelo
      tt <- terms(formula, data = data)
      # Crear matriz de diseño sin intercepto
      XTrain <- model.matrix(tt, data = data)[, -1, drop = FALSE]
      YTrain <- data$Clase
      # Entrenar modelo glmnet
      model <- glmnet(
        x = XTrain,
        y = YTrain,
        family = "multinomial",
        lambda = 0,
        type.multinomial = "ungrouped"
      )
      # Devolver modelo, términos y niveles de Clase
      list(
        glmnet_model = model,
        terms = tt,
        clase_levels = clase_levels
      )
    },
    predict.func = function(model, newdata) {
      # Extraer componentes del modelo
      glmnet_model <- model$glmnet_model
      tt <- model$terms
      clase_levels <- model$clase_levels
      # Crear matriz de diseño usando los términos del modelo de entrenamiento, sin intercepto
      tt_x <- delete.response(tt)
      XTest <- model.matrix(tt_x, data = newdata)
      # Eliminar intercepto
      XTest <- XTest[, -1, drop = FALSE]
      # Realizar predicción
      Pred_Class <- predict(
        object = glmnet_model,
        newx = XTest,
        type = "class"
      )
      # Convertir a factor con los niveles correctos
      Pred_Class <- factor(Pred_Class, levels = clase_levels)
      return(Pred_Class)
    },
    data.train = gen.data$Data$Train[[as.character(p)]],
    data.test = gen.data$Data$Test[[as.character(p)]]
  )
  
  # Almacenar resultados
  Models$MNL$Model[[as.character(p)]] <- result$models
  Models$MNL$Pred[[as.character(p)]] <- result$pred
  Models$MNL$MC[[as.character(p)]] <- result$MC
  Models$MNL$Accuracy[[as.character(p)]] <- result$accuracies
}


# Resultados --------------------------------------------------------------

# Construcción de la tabla de resultados incluyendo todos los modelos
Resultados <- bind_rows(
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$UGGM_QDA$Accuracy),
             Modelo = factor("UGGM-QDA")),
  
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$QDA$Accuracy),
             Modelo = factor("QDA")),
  
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$RF$Accuracy),
             Modelo = factor("RF")),
  
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$SVM$Accuracy),
             Modelo = factor("SVM")),
  
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$LDA$Accuracy),
             Modelo = factor("LDA")),
  
  data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
             N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
             Accuracy = unlist(Models$MNL$Accuracy),
             Modelo = factor("MLR"))

)

rownames(Resultados) <- NULL

# Ordenar los factores en la visualización
Resultados$N <- factor(Resultados$N, levels = c("50", "100", "500", "1000"))
Resultados$P <- factor(Resultados$P, levels = c("20", "25", "30", "35"))

# Generación de la gráfica
ggplot(data = Resultados, 
       mapping = aes(x = as.factor(N), 
                     y = Accuracy, 
                     color = Modelo, 
                     shape = Modelo)) +
  geom_point() +
  geom_line(aes(group = Modelo)) + 
  facet_grid(~P) + 
  theme_bw() + 
  labs(x = "Tamaño de muestra por clase", 
       y = "TCCG", 
       color = "Modelo", 
       shape = "Modelo",
       title = "Modelos de clasificación",
       subtitle = "Comparación de Modelos") +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))


