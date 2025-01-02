
###### Interval [1,4] PASS schedule --------------------------------------------
source("R script/Function/time dependent evaluation metrics.R")
load("R script/Simulation/Seed/DGM_seed.RData")

# Correctly specified model -----------------------------------------------
for (i in 1:200) {
  load(paste0("R script/Simulation/Models/Correctly specified models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("R script/Simulation/Datasets/traindata_", i+1, ".RData"))
  } else {
    load("R script/Simulation/Datasets/traindata_1.RData")
  }
  res <- td_eval(model = icjm, data = train.data, type = "corspe", 
                 t = 1, dt = 3, seed = seed[i])
  save(res, file=paste0("R script/Simulation/Evaluation/Correctly specified models/AccMet_", i, ".RData"))
}


# No covariate model ------------------------------------------------------------
for (i in 1:200) {
  load(paste0("R script/Simulation/Models/No covariate models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("R script/Simulation/Datasets/traindata_", i+1, ".RData"))
  } else {
    load("R script/Simulation/Datasets/traindata_1.RData")
  }
  try({res <- td_eval(model = icjm, data = train.data, type = "nocovar", 
                 t = 1, dt = 3, seed = seed[i])
  save(res, file=paste0("R script/Simulation/Evaluation/No covariate models/AccMet_", i, ".RData"))})
}

# Linear model ------------------------------------------------------------

for (i in 1:200) {
  load(paste0("R script/Simulation/Models/Linear models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("R script/Simulation/Datasets/traindata_", i+1, ".RData"))
  } else {
    load("R script/Simulation/Datasets/traindata_1.RData")
  }
  res <- td_eval(model = icjm, data = train.data, type = "linear", 
                 t = 1, dt = 3, seed = seed[i])
  save(res, file=paste0("R script/Simulation/Evaluation/Linear models/AccMet_", i, ".RData"))
}


###### Interval [1,4] Random schedule (0.3, 4) ------------------------------
source("R script/Function/time dependent evaluation metrics.R")
load("R script/Simulation/Seed/DGM_seed.RData")

# Correctly specified model -----------------------------------------------
for (i in 1:200) {
  try({
    load(paste0("R script/Simulation/Models/Correctly specified models/Joint model_", i, ".RData"))
    if (i != 200) {
      load(paste0("R script/Simulation/Datasets/random_3_40/traindata_", i+1, ".RData"))
    } else {
      load("R script/Simulation/Datasets/random_3_40/traindata_1.RData")
    }
    res <- td_eval(model = icjm, data = train.data, type = "corspe", 
                   t = 1, dt = 3, seed = seed[i])
    save(res, file=paste0("R script/Simulation/Evaluation/interval_1_4_random_3_40/Correctly specified models/AccMet_", i, ".RData"))
  })
}



###### Interval [1,4] Random schedule (0.3, 1) ------------------------------
source("R script/Function/time dependent evaluation metrics.R")
load("R script/Simulation/Seed/DGM_seed.RData")

# Correctly specified model -----------------------------------------------
for (i in 1:200) {
  try({
    load(paste0("R script/Simulation/Models/Correctly specified models/Joint model_", i, ".RData"))
    if (i != 200) {
      load(paste0("R script/Simulation/Datasets/random_3_10/traindata_", i+1, ".RData"))
    } else {
      load("R script/Simulation/Datasets/random_3_10/traindata_1.RData")
    }
    res <- td_eval(model = icjm, data = train.data, type = "corspe", 
                   t = 1, dt = 3, seed = seed[i])
    save(res, file=paste0("R script/Simulation/Evaluation/interval_1_4_random_3_10/Correctly specified models/AccMet_", i, ".RData"))
  })
}

###### Interval [1,4] Random schedule (1, 2) ------------------------------
source("R script/Function/time dependent evaluation metrics.R")
load("R script/Simulation/Seed/DGM_seed.RData")

# Correctly specified model -----------------------------------------------
for (i in 1:200) {
  try({
    load(paste0("R script/Simulation/Models/Correctly specified models/Joint model_", i, ".RData"))
    if (i != 200) {
      load(paste0("R script/Simulation/Datasets/random_10_20/traindata_", i+1, ".RData"))
    } else {
      load("R script/Simulation/Datasets/random_10_20/traindata_1.RData")
    }
    res <- td_eval(model = icjm, data = train.data, type = "corspe", 
                   t = 1, dt = 3, seed = seed[i])
    save(res, file=paste0("R script/Simulation/Evaluation/interval_1_4_random_10_20/Correctly specified models/AccMet_", i, ".RData"))
  })
}

