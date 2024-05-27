
###############
# Network QDA #
###############

library(igraph)
library(qgraph)
library(bootnet)
library(MASS)
library(dplyr)


# RAND --------------------------------------------------------------------

RAND <- function(alpha, p) {
  
  set.seed(123)
  
  network <- sample_pa(n = p, power = alpha, m = 1, directed = FALSE)
  
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


# Sample N(mu,Sigma) ------------------------------------------------------

Sample <- function(n, mu, Sigma, class.label) {
  
  mvrnorm(n = n, mu = mu, Sigma = Sigma) %>%
    data.frame() %>%
    mutate(Clase = as.factor(class.label))
  
}


# Funciones NetQDA --------------------------------------------------------

source("NetQDA.R")


# Simulaciones ------------------------------------------------------------

Gen.Data <- function(P, N) {
  
  Data <- list(Train = list(), Test = list())
  Network <- list()
  
  for (p in P) {
    
    RAND.1 <- RAND(alpha = 1, p = p)
    RAND.2 <- RAND(alpha = 1.1, p = p)
    RAND.3 <- RAND(alpha = 1.2, p = p)
    
    Network[[as.character(p)]] <- list(RAND.1$network, RAND.2$network, RAND.3$network)
    
    set.seed(123)
    data.train.p <- lapply(N, function(n) {
      muestra1.train <- Sample(n, RAND.1$m, RAND.1$S, "Clase1")
      muestra2.train <- Sample(n, RAND.2$m, RAND.2$S, "Clase2")
      muestra3.train <- Sample(n, RAND.3$m, RAND.3$S, "Clase3")
      bind_rows(muestra1.train, muestra2.train, muestra3.train)
    })
    
    names(data.train.p) <- as.character(N)
    Data$Train[[as.character(p)]] <- data.train.p
    
    set.seed(321)
    muestra1.test <- Sample(1000, RAND.1$m, RAND.1$S, "Clase1")
    muestra2.test <- Sample(1000, RAND.2$m, RAND.2$S, "Clase2")
    muestra3.test <- Sample(1000, RAND.3$m, RAND.3$S, "Clase3")
    Data$Test[[as.character(p)]] <- bind_rows(muestra1.test, muestra2.test, muestra3.test)
  }
  
  return(list(Network = Network, Data = Data))
  
}

P <- c(20, 25, 30, 35)
N <- c(50, 100, 500, 1000)

gen.data <- Gen.Data(P = P, N = N)


# Tuning ------------------------------------------------------------------

rho.tune <- list()

for (p in P) {
  
  rho.tune[[as.character(p)]] <- lapply(gen.data$Data$Train[[as.character(p)]], 
                                        function(Datos.Train) {
    tune.rho(formula = Clase~.,
             data = Datos.Train,
             rhos = seq(0, 1, by = 0.01),
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
  facet_grid(N~P) +
  theme_bw() +
  labs(x = "Rho", y = "Accuracy")


# NetQDA v.s QDA ----------------------------------------------------------

Aux <- function(model.func , predict.func, data.train, data.test, best.rho = NULL) {
  
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
    accuracies[[n]] <- sum(diag(mc))/sum(mc)
  }
  
  list(models = models, pred = pred, MC = mc, accuracies = accuracies)
}

Models <- list(NetQDA = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()),
               QDA = list(Model = list(), Pred = list(), MC = list(), Accuracy = list()))

for (p in P) {
  
  result <- Aux(model.func = NetQDA,
                predict.func = function(model, newdata) predict.NetQDA(object = model, newdata = newdata),
                data.train = gen.data$Data$Train[[as.character(p)]],
                data.test = gen.data$Data$Test[[as.character(p)]],
                best.rho = best.rho)
  
  Models$NetQDA$Model[[as.character(p)]] <- result$models
  Models$NetQDA$Pred[[as.character(p)]] <- result$pred
  Models$NetQDA$MC[[as.character(p)]] <- result$MC
  Models$NetQDA$Accuracy[[as.character(p)]] <- result$accuracies
  
}

for (p in P) {
  result <- Aux(model.func = qda,
                predict.func = function(model, newdata) predict(object = model, newdata = newdata)$class,
                data.train = gen.data$Data$Train[[as.character(p)]],
                data.test = gen.data$Data$Test[[as.character(p)]])
  
  Models$QDA$Model[[as.character(p)]] <- result$models
  Models$QDA$Pred[[as.character(p)]] <- result$pred
  Models$QDA$MC[[as.character(p)]] <- result$MC
  Models$QDA$Accuracy[[as.character(p)]] <- result$accuracies
}

Resultados <- bind_rows(data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
                                   N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
                                   Accuracy = unlist(Models$NetQDA$Accuracy),
                                   Modelo = factor("NetQDA")),
                        data.frame(P = rep(P, each = length(names(gen.data$Data$Train[[as.character(P[1])]]))),
                                   N = rep(names(gen.data$Data$Train[[as.character(P[1])]]), times = length(P)),
                                   Accuracy = unlist(Models$QDA$Accuracy),
                                   Modelo = factor("QDA")))

rownames(Resultados) <- NULL

Resultados$N <- factor(Resultados$N, levels=c("50", "100", "500", "1000"))
Resultados$P <- factor(Resultados$P, levels=c("20", "25", "30", "35"))

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
       y = "Accuracy", 
       color = "Modelo", 
       shape = "Modelo")


