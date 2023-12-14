
################
# Simulaciones #
################

library(dplyr)

library(MASS)

library(MVN)
library(mvnormtest)
library(mnt)
library(energy)

library(huge)

library(qgraph)
library(bootnet)

rm(list = ls(all.names = TRUE))
gc()

# AR(1) -------------------------------------------------------------------

S_AR1 <- function(rho, p){
  S <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      S[i,j] <- rho^(abs(i-j))
    }
  }
  return(S)
}

S1 <- S_AR1(0.5,50)
K1 <- solve(S1)

## Vector con valores de la media
m <- rep(x = 0, times = ncol(S1))

## Simulación de una muestra de tamaño n = 100
## con vector de medias "m" y matriz de varianzas
## y covarianzas "S1".
set.seed(789)
muestra1 <- mvrnorm(n = 1000, mu = m, Sigma = S1)

## Verificación Normal Multivariada.
mvn(data = muestra1, mvnTest = "mardia")$multivariateNormality
mvn(data = muestra1, mvnTest = "hz")$multivariateNormality
mvn(data = muestra1, mvnTest = "royston")$multivariateNormality

## mvnormtest (error numérico)
mshapiro.test(U = muestra1)

## mnt
test.BHEP(data = muestra1, 
          MC.rep = 100,
          alpha = 0.05)

## energy
mvnorm.test(x = muestra1, 
            R = 100)

## muestra1 (matriz) -> Muestra1 (DataFrame)
Muestra1 <- muestra1 %>% 
  data.frame()

## Modelo Gráfico Gaussiano con penalización lasso
MGG1 <- estimateNetwork(data = Muestra1,
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = 5)
plot(MGG1,
     edge.labels = FALSE,
     font = 2)

# MA(1) -------------------------------------------------------------------

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

S2 <- S_MA1(0.5,50)
K2 <- solve(S2)

## Vector con valores de la media
m <- rep(x = 0, times = ncol(S2))

## Simulación de una muestra de tamaño n = 100
## con vector de medias "m" y matriz de varianzas
## y covarianzas "S1".
set.seed(789)
muestra2 <- mvrnorm(n = 1000, mu = m, Sigma = S2)

## Verificación Normal Multivariada.
mvn(data = muestra2, mvnTest = "mardia")
mvn(data = muestra2, mvnTest = "hz")
mvn(data = muestra2, mvnTest = "royston")

## mvnormtest
mshapiro.test(U = muestra2)

## mnt
test.BHEP(data = muestra2, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]

## energy
mvnorm.test(x = muestra2, 
            R = 1000)

## muestra1 (matriz) -> Muestra1 (DataFrame)
Muestra2 <- muestra2 %>% 
  data.frame()

## Modelo Gráfico Gaussiano con penalización lasso
MGG2 <- estimateNetwork(data = Muestra2,
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = 1)
plot(MGG2,
     edge.labels = FALSE,
     font = 2)

Muestra2T <- Muestra2 %>% 
  exp()

## Verificación Normal Multivariada.
mvn(data = Muestra2T, mvnTest = "mardia")$multivariateNormality
mvn(data = Muestra2T, mvnTest = "hz")$multivariateNormality
mvn(data = Muestra2T, mvnTest = "royston")$multivariateNormality

MGG2T <- estimateNetwork(data = Muestra2T,
                         default = "EBICglasso",
                         corMethod = "cor",
                         tuning = 2)

plot(MGG2T,
     edge.labels = FALSE,
     font = 2)

Muestra2T_huge <- huge.npn(Muestra2T, npn.func = "truncation")

mvn(data = Muestra2T_huge, mvnTest = "mardia")$multivariateNormality
mvn(data = Muestra2T_huge, mvnTest = "hz")$multivariateNormality
mvn(data = Muestra2T_huge, mvnTest = "royston")$multivariateNormality

MGG2T_huge <- estimateNetwork(data = Muestra2T_huge,
                              default = "EBICglasso",
                              corMethod = "cor",
                              tuning = 5)

plot(MGG2T_huge,
     edge.labels = FALSE,
     font = 2)

# ECM ---------------------------------------------------------------------

S_ECM <- function(rho, p){
  S <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){
        S[i,j] <- 1
      } else{
        S[i,j] <- rho
      }
    }
  }
  return(S)
}

S3 <- S_ECM(0.5, 50)
K3 <- solve(S3)

## Vector con valores de la media
m <- rep(x = 0, times = ncol(S3))

## Simulación de una muestra de tamaño n = 100
## con vector de medias "m" y matriz de varianzas
## y covarianzas "S1".
set.seed(789)
muestra3 <- mvrnorm(n = 1000, mu = m, Sigma = S3)

## Verificación Normal Multivariada.
mvn(data = muestra3, mvnTest = "mardia")
mvn(data = muestra3, mvnTest = "hz")
mvn(data = muestra3, mvnTest = "royston")

## mvnormtest
mshapiro.test(U = muestra3)

## mnt
test.BHEP(data = muestra3, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]

## energy
mvnorm.test(x = muestra3, 
            R = 1000)

## muestra1 (matriz) -> Muestra1 (DataFrame)
Muestra3 <- muestra3 %>% 
  data.frame()

## Modelo Gráfico Gaussiano con penalización lasso
MGG3 <- estimateNetwork(data = Muestra3,
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = 0.5)
plot(MGG3,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)


# RAND --------------------------------------------------------------------

S_RAND <- function(p, R){
  set.seed(123)
  A <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      A[i,j] <- runif(1, min = -1, max = 1)
      }
    }
  diagonal <- c(NA)
  for (i in 1:p) {
    diagonal[i] <- sum(abs(A[i,]))-abs(A[i,i])
  }
  diag(A) <- R*diagonal
  AInv <- solve(A)
  B <- matrix(data = NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      B[i,j] <- AInv[i,j]/sqrt(AInv[i,i]*AInv[j,j])
    }
  }
  return(B)
}

S4 <- S_RAND(p=50, R=50)
K4 <- solve(S4)

## Vector con valores de la media
m <- rep(x = 0, times = ncol(S4))

## Simulación de una muestra de tamaño n = 100
## con vector de medias "m" y matriz de varianzas
## y covarianzas "S1".
set.seed(789)
muestra4 <- mvrnorm(n = 1000, mu = m, Sigma = S4)

## Verificación Normal Multivariada.
mvn(data = muestra4, mvnTest = "mardia")
mvn(data = muestra4, mvnTest = "hz")
mvn(data = muestra4, mvnTest = "royston")

## mvnormtest
mshapiro.test(U = muestra4)

## mnt
test.BHEP(data = muestra4, 
          MC.rep = 100,
          alpha = 0.05)[["Decision"]]

## energy
mvnorm.test(x = muestra4, 
            R = 1000)

## muestra1 (matriz) -> Muestra1 (DataFrame)
Muestra4 <- muestra4 %>% 
  data.frame()

## Modelo Gráfico Gaussiano con penalización lasso
MGG4 <- estimateNetwork(data = muestra4,
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = 0.5)
plot(MGG4,
     edge.labels = FALSE,
     font = 2)






###############
library(Matrix)
library(igraph)

sample <- function(l,s,sem){
  set.seed(sem)
  r <- c(0)
  for (i in 1:l) {
    r[i] <- runif(n = 1, min = 0, max = 1)
  }
  r <- sort(r)
  k <- 1
  bin <- c(NULL)
  for (i in 1:l) {
    while (r[i]>=s[k]) {
      k <- k+1
    }
    bin[i] <- k
  }
  return(bin)
}

Chung_Lu <- function(w){
  n <- length(w)
  wsum <- c(NULL)
  for (i in 1:n) {
    if(i == 1){
      wsum[i] <- 0 + w[i]
    }else{
      wsum[i] <- wsum[i-1] + w[i]
    }
  }
  wsum <- wsum/wsum[n]
  l <- ceiling(0.5*(sum(w)+mean(w)^2))
  I <- sample(l, wsum, sem = 123)
  J <- sample(l, wsum, sem = 789)
  M <- get.adjacency(graph.edgelist(matrix(c(I,J), 
                                           ncol = 2), 
                                    directed = FALSE))
  A <- as.matrix(M)
  diag(A) <- 0
  n <- dim(M)[1]
  D <- matrix(0, ncol = n,
                 nrow = n)
  for (i in 1:n) {
    D[i,i] <- sum(M[i,])
  }
  #M <- +as.matrix(sparseMatrix(I,J))
  #M <- get.adjacency(graph.edgelist(as.matrix(sparseMatrix(I,J))))
  #M <- +as.matrix(sparseMatrix(I,J))
  #M <- (M + t(M))
  #diag(M) <- rep(1,dim(M)[1])
  return(D-A)
}

PLRw <- function(beta,w,alfa){
  n <- length(w)
  c <- mean(w)*((beta-2)/(beta-1))
  wM <- n^alfa # max(w)
  i0 <- n*(c/wM)^(beta-1)
  wPLR <- c(NULL)
  for (i in 1:n) {
    wPLR[i] <- c*((i+i0-1)/n)^(-1/(beta-1))
  }
  return(wPLR)
}

q <- PLRw(beta = 2.5,
          alfa = 0.45,
          w = seq(from = 1, to = 20, by = 0.1))

M <- Chung_Lu(q)

K <- diag(dim(M)[1]) + M
S <- solve(K)

## Vector con valores de la media
m <- rep(x = 0, times = ncol(S))

## Simulación de una muestra de tamaño n = 100
## con vector de medias "m" y matriz de varianzas
## y covarianzas "S1".
set.seed(789)
muestra4 <- mvrnorm(n = 1000, mu = m, Sigma = S)

## muestra1 (matriz) -> Muestra1 (DataFrame)
Muestra4 <- muestra4 %>% 
  as.data.frame()

## Modelo Gráfico Gaussiano con penalización lasso
MGG4 <- estimateNetwork(data = muestra4,
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = 1.5)
plot(MGG4,
     edge.labels = FALSE,
     font = 2)

