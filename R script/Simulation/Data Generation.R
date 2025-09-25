library(PCaASSim)

set.seed(2024)
seed <- sample(1:1e5, 200)
save(seed, file="R script/Simulation/Seed/DGM_seed.RData")

for (i in 1:200) {
  sim <- icjmsim(n = 1000, seed = seed[i])
  train.data <- sim$dat[sim$dat$CISNET_ID %in% 1:300,]
  train.data.id <- sim$dat.id[sim$dat.id$CISNET_ID %in% 1:300,]
  save(train.data,
       file = paste0("Data/Simulation Datasets/traindata_",
                     i,
                     ".RData"))
  save(train.data.id,
       file = paste0("Data/Simulation Datasets/traindata_id_",
                     i,
                     ".RData"))
}


# Generating random biopsy intervals [0.3,4] --------------------------------------
set.seed(2025)
seed <- sample(1:1e5, 200)

for (i in 1:200) {
  load(paste0("Data/Simulation Datasets/traindata_id_", i, ".RData"))
  load(paste0("Data/Simulation Datasets/traindata_", i, ".RData"))
  # his schedules
  set.seed(seed[i])
  n <- nrow(train.data.id)
  cmp1 <- rep(NA, n)
  cmp2 <- rep(NA, n)
  status <- rep(NA, n)
  for (j in 1:n) {
    intervals <- c(0, runif(50, 0.3, 4))
    schedules <- cumsum(intervals)
    t.det <- schedules[schedules>=train.data.id$time.prg[j]][1]
    t.trt <- train.data.id$time.trt[j]
    t.cen <- train.data.id$time.cen[j]
    # status.cmp
    status[j] <- which.min(c(t.cen, t.det, t.trt)) - 1
    # cmp1
    cmp1[j] <- max(schedules[schedules<min(t.cen, t.det, t.trt)])
    # cmp2
    cmp2[j] <- min(t.cen, t.det, t.trt)
  }
  train.data.id$time.cmp1 <- cmp1
  train.data.id$time.cmp2 <- cmp2
  train.data.id$status.cmp <- status
  save(train.data.id,
       file = paste0("Data/Simulation Datasets/random_3_40/traindata_id_",
                     i,
                     ".RData"))
  
  ### long dataset
  nfreq <- table(train.data$CISNET_ID)
  cmp1.long <- rep(cmp1, nfreq)
  cmp2.long <- rep(cmp2, nfreq)
  status.long <- rep(status, nfreq)
  train.data$time.cmp1 <- cmp1.long
  train.data$time.cmp2 <- cmp2.long
  train.data$status.cmp <- status.long
  save(train.data,
       file = paste0("Data/Simulation Datasets/random_3_40/traindata_",
                     i,
                     ".RData"))
}

# Generating random biopsy intervals [0.3,1] --------------------------------------
set.seed(2025)
seed <- sample(1:1e5, 200)

for (i in 1:200) {
  load(paste0("Data/Simulation Datasets/traindata_id_", i, ".RData"))
  load(paste0("Data/Simulation Datasets/traindata_", i, ".RData"))
  # his schedules
  set.seed(seed[i])
  n <- nrow(train.data.id)
  cmp1 <- rep(NA, n)
  cmp2 <- rep(NA, n)
  status <- rep(NA, n)
  for (j in 1:n) {
    intervals <- c(0, runif(50, 0.3, 1))
    schedules <- cumsum(intervals)
    t.det <- schedules[schedules>=train.data.id$time.prg[j]][1]
    t.trt <- train.data.id$time.trt[j]
    t.cen <- train.data.id$time.cen[j]
    # status.cmp
    status[j] <- which.min(c(t.cen, t.det, t.trt)) - 1
    # cmp1
    cmp1[j] <- max(schedules[schedules<min(t.cen, t.det, t.trt)])
    # cmp2
    cmp2[j] <- min(t.cen, t.det, t.trt)
  }
  train.data.id$time.cmp1 <- cmp1
  train.data.id$time.cmp2 <- cmp2
  train.data.id$status.cmp <- status
  save(train.data.id,
       file = paste0("Data/Simulation Datasets/random_3_10/traindata_id_",
                     i,
                     ".RData"))
  
  ### long dataset
  nfreq <- table(train.data$CISNET_ID)
  cmp1.long <- rep(cmp1, nfreq)
  cmp2.long <- rep(cmp2, nfreq)
  status.long <- rep(status, nfreq)
  train.data$time.cmp1 <- cmp1.long
  train.data$time.cmp2 <- cmp2.long
  train.data$status.cmp <- status.long
  save(train.data,
       file = paste0("Data/Simulation Datasets/random_3_10/traindata_",
                     i,
                     ".RData"))
}


# Generating random biopsy intervals [1,2] --------------------------------------
set.seed(2025)
seed <- sample(1:1e5, 200)

for (i in 1:200) {
  load(paste0("Data/Simulation Datasets/traindata_id_", i, ".RData"))
  load(paste0("Data/Simulation Datasets/traindata_", i, ".RData"))
  # his schedules
  set.seed(seed[i])
  n <- nrow(train.data.id)
  cmp1 <- rep(NA, n)
  cmp2 <- rep(NA, n)
  status <- rep(NA, n)
  for (j in 1:n) {
    intervals <- c(0, runif(50, 1, 2))
    schedules <- cumsum(intervals)
    t.det <- schedules[schedules>=train.data.id$time.prg[j]][1]
    t.trt <- train.data.id$time.trt[j]
    t.cen <- train.data.id$time.cen[j]
    # status.cmp
    status[j] <- which.min(c(t.cen, t.det, t.trt)) - 1
    # cmp1
    cmp1[j] <- max(schedules[schedules<min(t.cen, t.det, t.trt)])
    # cmp2
    cmp2[j] <- min(t.cen, t.det, t.trt)
  }
  train.data.id$time.cmp1 <- cmp1
  train.data.id$time.cmp2 <- cmp2
  train.data.id$status.cmp <- status
  save(train.data.id,
       file = paste0("Data/Simulation Datasets/random_10_20/traindata_id_",
                     i,
                     ".RData"))
  
  ### long dataset
  nfreq <- table(train.data$CISNET_ID)
  cmp1.long <- rep(cmp1, nfreq)
  cmp2.long <- rep(cmp2, nfreq)
  status.long <- rep(status, nfreq)
  train.data$time.cmp1 <- cmp1.long
  train.data$time.cmp2 <- cmp2.long
  train.data$status.cmp <- status.long
  save(train.data,
       file = paste0("Data/Simulation Datasets/random_10_20/traindata_",
                     i,
                     ".RData"))
}
