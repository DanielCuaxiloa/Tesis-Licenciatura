########################
# Funciones Auxiliares # 
########################

rm(list = ls(all.names = TRUE))
gc()


# Funciones auxiliares ----------------------------------------------------

ConjuntoEvaluacion <- function(x) {
  
  Train <- training(x)
  Test <- testing(x)
  
  ConjuntoEvaluacion <- list(Train = Train,
                             Test = Test) 
  
  return(ConjuntoEvaluacion)
  
}

ErroresClasificacion <- function(MC.Train, MC.Test){
  
  # Train
  
  ## Clase 1
  TrainClase1 <- MC.Train[1,1]/sum(MC.Train[1,])
  
  ## Clase 2
  TrainClase2 <- MC.Train[2,2]/sum(MC.Train[2,])
  
  ## Clase 3
  TrainClase3 <- MC.Train[3,3]/sum(MC.Train[3,])
  
  ## TCC (Accuracy)
  TrainGlobal <- sum(diag(MC.Train))/sum(MC.Train)
  
  # Test
  
  ## Clase 1
  TestClase1 <- MC.Test[1,1]/sum(MC.Test[1,])
  
  ## Clase 2
  TestClase2 <- MC.Test[2,2]/sum(MC.Test[2,])
  
  ## Clase 3
  TestClase3 <- MC.Test[3,3]/sum(MC.Test[3,])
  
  ## TCC (Accuracy)
  TestGlobal <- sum(diag(MC.Test))/sum(MC.Test)
  
  Errores <- data.frame(TrainClase1, TrainClase2, TrainClase3, TrainGlobal,
                        TestClase1,  TestClase2,  TestClase3,  TestGlobal)
  
  return(Errores)
  
}

Evaluacion <- function(Metodo, workers){

  plan(strategy = multisession, 
       workers = workers)
  
  set.seed(1234)
  
  Individual <- map(.x = Folds$splits,
                    .f = ~ConjuntoEvaluacion(.x)) %>% 
    transpose() %>% 
    future_pmap(.f = get(Metodo),
                .options = furrr_options(seed = TRUE),
                .progress = TRUE) %>% 
    transpose() %>% 
    pmap(.f = ~ErroresClasificacion(.x,.y)) %>% 
    ldply(data.frame) %>% 
    mutate_if(is.numeric, ~.*100) %>% 
    mutate_if(is.numeric, round, 3) %>% 
    add_column(Metodo = Metodo,
               .before = "TrainClase1")
  
  Global <- Individual %>% 
    mutate(ID1 = factor(Folds$id)) %>% 
    group_by(ID1) %>% 
    summarise(TrainClase1 = mean(TrainClase1),
              TrainClase2 = mean(TrainClase2),
              TrainClase3 = mean(TrainClase3),
              TrainGlobal = mean(TrainGlobal),
              TestClase1 = mean(TestClase1),
              TestClase2 = mean(TestClase2),
              TestClase3 = mean(TestClase3),
              TestGlobal = mean(TestGlobal))
  
  return(list(Individual = Individual,
              Global = Global))
  
}

