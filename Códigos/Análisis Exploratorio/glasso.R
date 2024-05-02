
##########
# glasso #
##########

library(glasso)
library(bootnet)
library(NetDA)

library(MASS)
library(dplyr)


# Estructuras -------------------------------------------------------------

AR1 <- function(rho, p){
  
  S <- matrix(data = NA, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      S[i,j] <- rho^(abs(i-j))
    }
  }
  
  return(S)
  
}

S_MA1 <- function(rho, p){
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


# Simulaciones ------------------------------------------------------------

S1 <- AR1(rho = 0.25, p = 10)
K1 <- solve(S1)
m1 <- rep(x = 0, times = ncol(S1))

S2 <- S_MA1(rho = 0.5, p = 10)
K2 <- solve(S2)
m2 <- rep(x = 3, times = ncol(S2))

set.seed(789)

muestra1 <- mvrnorm(n = 1000, mu = m1, Sigma = S1) %>% 
  data.frame() %>% 
  mutate(Clase = 1)

muestra2 <- mvrnorm(n = 1000, mu = m2, Sigma = S2) %>% 
  data.frame() %>% 
  mutate(Clase = 2)


# glasso ------------------------------------------------------------------

Network1 <- glasso(s = K1,
                   rho = 0.1)

Network1$wi


# glasso LDA --------------------------------------------------------------

Datos <- bind_rows(muestra1, muestra2) 

X <- Datos %>% 
  dplyr::select(-Clase) %>% 
  as.matrix()

Y <- Datos %>%
  dplyr::select(Clase) %>% 
  as.matrix()

rho <- 0.01
class <- max(Y)
n <- length(Y)
p <- dim(X)[2]
pi <- NULL
x <- NULL
meanx <- NULL

## Probabilidades a priori
## vector de medias
for (i in 1:class) {
  pi <- c(pi, length(which(Y == i))/n)
  x[[i]] <- X[which(Y == i), ]
  meanx[[i]] <- as.vector(colMeans(x[[i]]))
}

## Matriz de precisión para cada clase (QDA)
Q_GLASSO <- NULL
for (i in 1:class) {
  s <- as.matrix(var(x[[i]]))
  Theta_glasso <- glasso(s, rho = rho)$wi
  Q_GLASSO[[i]] <- Theta_glasso
}

## Matriz de precisión asumiendo igualda (LDA)
ls <- as.matrix(var(X))
L_GLASSO <- glasso(ls, rho = rho)$wi

## Predicciones
predClass <- NULL
X_test <- X
n1 <- dim(X_test)[1]

## LDA
{
  EstPrecision <- L_GLASSO
  Cmatrix <- NULL
  for (i in 1:class) {
    Cmatrix <- cbind(Cmatrix, log(pi[i]) - 0.5 * as.numeric(t(meanx[[i]]) %*% 
                                                             L_GLASSO %*% meanx[[i]]) + as.vector(meanx[[i]] %*% 
                                                                                                    L_GLASSO %*% t(X_test)))
  }
  for (j in 1:n1) {
    predClass <- c(predClass, which(Cmatrix[j, ] == max(Cmatrix[j,])))
  }
}

## QDA
{
  EstPrecision <- Q_GLASSO
  Cmatrix <- NULL
  for (i in 1:class) {
    Cmatrix <- cbind(Cmatrix, log(pi[i]) - 0.5 * log(det(Q_GLASSO[[i]])) - 
                      0.5 * diag(t(as.matrix(t(X_test) - meanx[[i]])) %*% 
                                   Q_GLASSO[[i]] %*% (as.matrix(t(X_test) - meanx[[i]]))))
  }
  for (j in 1:n1) {
    predClass = c(predClass, which(Cmatrix[j, ] == max(Cmatrix[j, 
    ])))
  }
}

Modelo1 <- NetDA(X = X,
                 Y = Y,
                 method = 1,
                 X_test = X)

table(Y,Modelo1$yhat)

Modelo2 <- NetDA(X = X,
                 Y = Y,
                 method = 2,
                 X_test = X)

table(Y,Modelo2$yhat)

# EBICglasso --------------------------------------------------------------

Network1 <- estimateNetwork(data = dplyr::select(muestra1, -Clase),
                            default = "EBICglasso",
                            corMethod = "cor_auto",
                            tuning = 0.5)

Network2 <- estimateNetwork(data = dplyr::select(muestra2, -Clase),
                            default = "EBICglasso",
                            corMethod = "cor_auto",
                            tuning = 0.5)

plot(Network1,
     layout = "circle",
     edge.labels = FALSE,
     label.prop = 1,
     font = 2)

plot(Network2,
     layout = "circle",
     edge.labels = FALSE,
     label.prop = 1,
     font = 2)




