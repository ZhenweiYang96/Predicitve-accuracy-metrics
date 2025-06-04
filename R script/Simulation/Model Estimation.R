# load parallel packages -------------------------------------
library(future)
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
load("R script/Simulation/Seed/DGM_seed.RData")


# correctly-specified models -------------------------------------
source("R script/Function/ICJM_correctlyspecify.R")
plan(
  list(
    tweak(multisession, workers = 5),
    tweak(multisession, workers = 3)
  )
)

t1 <- Sys.time()
out <- lapply(1:200, function(datanum) {
  future({
    
    # load data file
    load(paste0("Data/Simulation Datasets/traindata_",
                datanum, ".RData"))
    # load inits
    load(paste0("Output/Simulation Models/Fitted mixed model/MM_inits_",
                datanum, ".RData"))
    
    # start modelling
    icjm <- icjm_cor(data = train.data, n.adapt = 3000, 
                     n.burnin = 3000, n.iter = 10000, 
                     GLMM_fit = mm_inits,
                     seed = seed[datanum])
    save(icjm,
         file = paste0("Output/Simulation Models/Correctly specified models/Joint model_",
                       datanum, ".RData"))
    # 
    plan(
      list(
        tweak(multisession, workers = 5),
        tweak(multisession, workers = 3)
      )
    )
    return(paste0(datanum,"done!"))
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()

t2 - t1



# no baseline covariate models -------------------------------------
source("R script/Function/ICJM_nocovar.R")
plan(
  list(
    tweak(multisession, workers = 5),
    tweak(multisession, workers = 3)
  )
)

t1 <- Sys.time()
out <- lapply(1:200, function(datanum) {
  future({
    
    # load data file
    load(paste0("Data/Simulation Datasets/traindata_",
                datanum, ".RData"))
    # load inits
    load(paste0("Output/Simulation Models/Fitted mixed model/MM_inits_",
                datanum, ".RData"))
    
    # start modelling
    icjm <- icjm_nocorv(data = train.data, n.adapt = 3000, 
                     n.burnin = 3000, n.iter = 10000, 
                     GLMM_fit = mm_inits,
                     seed = seed[datanum])
    save(icjm,
         file = paste0("Output/Simulation Models/No covariate models/Joint model_",
                       datanum, ".RData"))
    # 
    plan(
      list(
        tweak(multisession, workers = 5),
        tweak(multisession, workers = 3)
      )
    )
    return(paste0(datanum,"done!"))
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()

t2 - t1


# linear models -------------------------------------
source("Function/ICJM_linear.R")
plan(
  list(
    tweak(multisession, workers = 5),
    tweak(multisession, workers = 3)
  )
)

t1 <- Sys.time()
out <- lapply(1:100, function(datanum) {
  future({
    
    # load data file
    load(paste0("Data/Simulation Datasets/traindata_",
                datanum, ".RData"))
    
    # start modelling
    icjm <- icjm_linear(data = train.data, n.adapt = 3000, 
                        n.burnin = 3000, n.iter = 10000, 
                        seed = seed[datanum])
    save(icjm,
         file = paste0("Output/Simulation Models/Linear models/Joint model_",
                       datanum, ".RData"))
    # 
    plan(
      list(
        tweak(multisession, workers = 5),
        tweak(multisession, workers = 3)
      )
    )
    # return(paste0(datanum,"done!"))
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()

t2 - t1
