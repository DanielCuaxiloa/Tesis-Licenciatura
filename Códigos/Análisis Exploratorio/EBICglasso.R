
######################
# Modelos EBICglasso #      
######################


library(purrr)

# Modelos -----------------------------------------------------------------

ebicMin <- function(Grupo, alfa) {
  
  EN <- estimateNetwork(data = Grupo[,Variables],
                        default = "EBICglasso",
                        corMethod = "cor",
                        tuning = alfa)
  
  ebic_min <- EN[["results"]][["ebic"]] %>% min()
  
  return(list(EBIC = ebic_min,
              Alfa = alfa))
  
}

alfa <- seq(from = 0, to = 5, by = 0.1)

M1 <- map(.x = alfa, 
         .f = ~ebicMin(Grupo = Grupo1T, alfa = .x)) %>% 
  transpose()

M2 <- map(.x = alfa, 
          .f = ~ebicMin(Grupo = Grupo2T, alfa = .x)) %>% 
  transpose()

M3 <- map(.x = alfa, 
          .f = ~ebicMin(Grupo = Grupo3T, alfa = .x)) %>% 
  transpose()

Alfa_M1T <- unlist(M1[["Alfa"]])[M1[["EBIC"]] %>% which.min()]
Alfa_M2T <- unlist(M2[["Alfa"]])[M2[["EBIC"]] %>% which.min()]
Alfa_M3T <- unlist(M3[["Alfa"]])[M3[["EBIC"]] %>% which.min()]


# Grupo 1 -----------------------------------------------------------------

Modelo1T <- estimateNetwork(data = Grupo1T[,Variables],
                            default = "EBICglasso",
                            corMethod = "cor",
                            tuning = Alfa_M1T)
plot(Modelo1T,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

Modelo1 <- estimateNetwork(data = Grupo1[,Variables],
                           default = "EBICglasso",
                           corMethod = "cor",
                           tuning = Alfa_M1T)
plot(Modelo1,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

# Grupo 2 -----------------------------------------------------------------

Modelo2T <- estimateNetwork(data = Grupo2T[,Variables],
                            default = "EBICglasso",
                            corMethod = "cor",
                            tuning = Alfa_M2T)

plot(Modelo2T,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)

Modelo2 <- estimateNetwork(data = Grupo2[,Variables],
                           default = "EBICglasso",
                           corMethod = "cor",
                           tuning = Alfa_M2T)
plot(Modelo2,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)


# Grupo 3 -----------------------------------------------------------------

Modelo3T <- estimateNetwork(data = Grupo3T[,Variables],
                            default = "EBICglasso",
                            corMethod = "cor",
                            tuning = Alfa_M3T)

plot(Modelo3T,
     layout = "circle",
     edge.labels = FALSE,
     font = 0.5)

Modelo3 <- estimateNetwork(data = Grupo3[,Variables],
                           default = "EBICglasso",
                           corMethod = "cor",
                           tuning = Alfa_M3T)
plot(Modelo3,
     layout = "circle",
     edge.labels = FALSE,
     font = 2)



