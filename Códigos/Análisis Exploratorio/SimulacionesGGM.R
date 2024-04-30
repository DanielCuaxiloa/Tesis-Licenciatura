

library(MASS)
library(dplyr)
library(NetDA)

AR1 <- function(rho, p){
  
  S <- matrix(data = NA, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      S[i,j] <- rho^(abs(i-j))
    }
  }
  
  return(S)

  }

S1 <- AR1(rho = 0.25, p = 20)
S2 <- AR1(rho = 0.50, p = 20)
S3 <- AR1(rho = 0.75, p = 20)

K1 <- solve(S1)
K2 <- solve(S2)
K3 <- solve(S3)

m1 <- rep(x = 5, times = ncol(S1))
m2 <- rep(x = 7, times = ncol(S1))
m3 <- rep(x = 10, times = ncol(S1))

set.seed(789)

muestra1 <- mvrnorm(n = 1000, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase1"))

muestra2 <- mvrnorm(n = 1000, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase2"))

muestra3 <- mvrnorm(n = 1000, mu = m3, Sigma = S3) %>% 
  data.frame() %>% 
  mutate(Clase = as.factor("Clase3"))

Datos <- bind_rows(muestra1, muestra2, muestra3)


# Modelo LDA --------------------------------------------------------------

Modelo1 <- lda(formula = Clase ~.,
               data = Datos)

Pred.Modelo1 <- predict(object = Modelo1,
                        newdata = Datos)$class

table(Datos$Clase, Pred.Modelo1)


# NetLDA ------------------------------------------------------------------

X <- Datos %>% 
  dplyr::select(-Clase) %>% 
  as.matrix()

Y <- Datos %>% 
  mutate(Clase = case_when(
    Clase == "Clase1" ~ 1,
    Clase == "Clase2" ~ 2,
    Clase == "Clase3" ~ 3)) %>% 
  dplyr::select(Clase) %>% 
  as.matrix()

Modelo1.1 <- NetDA(X = X,
                   Y = Y,
                   method = 1,
                   X_test = X)

table(Y,Modelo1.1$yhat)


# Modelo QDA --------------------------------------------------------------

Modelo2 <- qda(formula = Clase ~.,
               data = Datos)

Pred.Modelo2 <- predict(object = Modelo2,
                        newdata = Datos)$class

table(Datos$Clase, Pred.Modelo2)


# NetQDA ------------------------------------------------------------------

Modelo2.1 <- NetDA(X = X,
                   Y = Y,
                   method = 2,
                   X_test = X)

table(Y,Modelo2.1$yhat) 
