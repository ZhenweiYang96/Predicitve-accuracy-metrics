#################################################
### True Brier scores and AUC
#################################################

# -----------------------------------------------------------------
# The true AUC and Brier score is calculated in absence of censoring
# The predicted risk at t+dt is predicted from the correctly-specified model
# -----------------------------------------------------------------

library(tidyverse)

for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Correctly specified models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    group[j] <- ifelse(df$time.prg[j] < df$time.trt[j] & 
                         df$time.prg[j] >= df$t[j] & 
                         df$time.prg[j] < df$tdt[j], "case", "control")
  }
  df$group <- group
  
  ## calculate AUC
  df_case <- df[df$group == "case",]
  df_control <- df[df$group == "control",]
  
  auc_df <- expand.grid(df_case$CISNET_ID, df_control$CISNET_ID)
  colnames(auc_df) <- c("ID_case", "ID_control")
  auc_df$auc_comp <- NA
  for (j in 1:nrow(auc_df)) {
    auc_df$auc_comp[j] <- ifelse(df_case$risk_tdt[df_case$CISNET_ID == auc_df$ID_case[j]] >
                                   df_control$risk_tdt[df_control$CISNET_ID == auc_df$ID_control[j]], 1, 0)
  }
  true_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2
  })
  true_bs <- mean(bs_df$sq)
  
  true_res <- list(auc_df = auc_df,
                   true_auc = true_auc,
                   bs_df = bs_df,
                   true_bs = true_bs)
  save(true_res, file = paste0("Output/Simulation Evaluation/True models/AccMet_true_",
                               i, ".RData"))
}


# True evaluation of no covariate model -----------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/No covariate models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    group[j] <- ifelse(df$time.prg[j] < df$time.trt[j] & 
                         df$time.prg[j] >= df$t[j] & 
                         df$time.prg[j] < df$tdt[j], "case", "control")
  }
  df$group <- group
  
  ## calculate AUC
  df_case <- df[df$group == "case",]
  df_control <- df[df$group == "control",]
  
  auc_df <- expand.grid(df_case$CISNET_ID, df_control$CISNET_ID)
  colnames(auc_df) <- c("ID_case", "ID_control")
  auc_df$auc_comp <- NA
  for (j in 1:nrow(auc_df)) {
    auc_df$auc_comp[j] <- ifelse(df_case$risk_tdt[df_case$CISNET_ID == auc_df$ID_case[j]] >
                                   df_control$risk_tdt[df_control$CISNET_ID == auc_df$ID_control[j]], 1, 0)
  }
  true_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2
  })
  true_bs <- mean(bs_df$sq)
  
  true_res <- list(auc_df = auc_df,
                   true_auc = true_auc,
                   bs_df = bs_df,
                   true_bs = true_bs)
  save(true_res, file = paste0("Output/Simulation Evaluation/True no cov/AccMet_true_",
                               i, ".RData"))
}


# True evaluation of linear model -----------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Linear models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    group[j] <- ifelse(df$time.prg[j] < df$time.trt[j] & 
                         df$time.prg[j] >= df$t[j] & 
                         df$time.prg[j] < df$tdt[j], "case", "control")
  }
  df$group <- group
  
  ## calculate AUC
  df_case <- df[df$group == "case",]
  df_control <- df[df$group == "control",]
  
  auc_df <- expand.grid(df_case$CISNET_ID, df_control$CISNET_ID)
  colnames(auc_df) <- c("ID_case", "ID_control")
  auc_df$auc_comp <- NA
  for (j in 1:nrow(auc_df)) {
    auc_df$auc_comp[j] <- ifelse(df_case$risk_tdt[df_case$CISNET_ID == auc_df$ID_case[j]] >
                                   df_control$risk_tdt[df_control$CISNET_ID == auc_df$ID_control[j]], 1, 0)
  }
  true_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2
  })
  true_bs <- mean(bs_df$sq)
  
  true_res <- list(auc_df = auc_df,
                   true_auc = true_auc,
                   bs_df = bs_df,
                   true_bs = true_bs)
  save(true_res, file = paste0("Output/Simulation Evaluation/True linear/AccMet_true_",
                               i, ".RData"))
}


#################################################
### True EPCE
#################################################
library(tidyverse)
source("R script/Function/time dependent true EPCE.R")
load("R script/Simulation/Seed/DGM_seed.RData")

# Correct model
for (i in 1:200) {
  load(paste0("Output/Simulation Models/Correctly specified models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("Data/Simulation Datasets/traindata_", i+1, ".RData"))
  } else {
    load("Data/Simulation Datasets/traindata_1.RData")
  }
  true_res <- true_epce(model = icjm, data = train.data, type = "corspe", 
                   t = 1, dt = 3, seed = seed[i])
  save(true_res, file=paste0("Output/Simulation Evaluation/True models/EPCE_true_", i, ".RData"))
}

# Linear model
for (i in 1:200) {
  load(paste0("Output/Simulation Models/Linear models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("Data/Simulation Datasets/traindata_", i+1, ".RData"))
  } else {
    load("Data/Simulation Datasets/traindata_1.RData")
  }
  true_res <- true_epce(model = icjm, data = train.data, type = "linear", 
                        t = 1, dt = 3, seed = seed[i])
  save(true_res, file=paste0("Output/Simulation Evaluation/True linear/EPCE_true_", i, ".RData"))
}

# No baseline covariate model
for (i in 1:200) {
  load(paste0("Output/Simulation Models/No covariate models/Joint model_", i, ".RData"))
  if (i != 200) {
    load(paste0("Data/Simulation Datasets/traindata_", i+1, ".RData"))
  } else {
    load("Data/Simulation Datasets/traindata_1.RData")
  }
  true_res <- true_epce(model = icjm, data = train.data, type = "nocovar", 
                        t = 1, dt = 3, seed = seed[i])
  save(true_res, file=paste0("Output/Simulation Evaluation/True no cov/EPCE_true_", i, ".RData"))
}

