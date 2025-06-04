#################################################
### Observed Brier scores and AUC
#################################################

# -----------------------------------------------------------------
# The observed AUC and Brier score is calculated using the observed progression time
# The predicted risk at t+dt is predicted from the corresponding model
# -----------------------------------------------------------------

library(tidyverse)

# Observed evaluation of correctly-specified model ----------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Correctly specified models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed models/AccMet_obs_",
                              i, ".RData"))
}


# Observed evaluation of no covariate model -----------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/No covariate models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed no cov/AccMet_obs_",
                              i, ".RData"))
}


# Observed evaluation of linear model -----------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Linear models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed linear/AccMet_obs_",
                              i, ".RData"))
}


# Observed evaluation of random schedule [0.3, 1] -------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_3_10/Correctly specified models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_3_10/AccMet_obs_",
                              i, ".RData"))
}


# Observed evaluation of random schedule [1, 2] ---------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_10_20/Correctly specified models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_10_20/AccMet_obs_",
                              i, ".RData"))
}

# Observed evaluation of random schedule [0.3, 4] -------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_3_40/Correctly specified models/AccMet_", 
              i, ".RData"))
  df <- res$data.id
  group <- rep(NA, nrow(df))
  for (j in 1:nrow(df)) {
    if (df$status.cmp[j] == 1 & 
        df$time.cmp2[j] >= df$t[j] & 
        df$time.cmp2[j] < df$tdt[j]) {
      group[j] <- "case"
    } else if (df$time.cmp2[j] > df$tdt[j]) {
      group[j] <- "control"
    } else {
      group[j] <- "excluded"
    }
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
  obs_auc <- mean(auc_df$auc_comp)
  
  ## calculate Brier scores
  bs_df <- df[,c("CISNET_ID", "risk_tdt", "group")]
  bs_df$sq <- sapply(1:nrow(bs_df), function(x) {
    ifelse(bs_df$group[x] == "excluded", 0, (bs_df$risk_tdt[x] - (bs_df$group[x] == "case"))^2)
  })
  obs_bs <- mean(bs_df$sq)
  
  obs_res <- list(auc_df = auc_df,
                  obs_auc = obs_auc,
                  bs_df = bs_df,
                  obs_bs = obs_bs)
  save(obs_res, file = paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_3_40/AccMet_obs_",
                              i, ".RData"))
}

