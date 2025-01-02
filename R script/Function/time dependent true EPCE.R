##################################################
### Function for evaluation accuracy metrics 
### (time-dependent AUC, Brier score, and EPCE)
##################################################

# library(PCaASSim)
library(survival)
library(future)
library(parallel)
# load("Data/pass.RData")
# load("Model/ICJM1.RData")
# load("Simulation/Models/No covariate models/Joint model_1.RData")
# load("Simulation/Models/Linear models/Joint model_1.RData")
# model <- icjm#ICJM1
# data <- pass
# data <- data[data$CISNET_ID < 50,]
# eidVar = "status.cmp"
# t=1
# dt=3
# pred_niter=200
# seed = 2024
# c.step=50
# qp = 15
# type = "linear"
# type = "nocovar"
# source("Function/dynpred_cor.R")
# source("Function/function.R")

true_epce <- function(model, data, 
                      type = "corspe",
                      qp = 15,
                      t, dt, eidVar = "status.cmp", pred_niter=200,
                      seed = 2024, c.step=50 ) { # , infty = 0, idVar = "CISNET_ID", timeVar = "TimeSince_Dx",e1tVar = "time.cmp1", e2tVar = "time.cmp2",
  if (!(type %in% c('corspe', 'nocovar','linear'))) {
    stop("there are three models to be chosen from, you should specify 'corspe', 'nocovar', or 'linear'.")
  }
  source(paste0("Function/dynpred_", type, ".R"))
  
  # variable names
  idVar <- model$model_info$var_names$id
  longitime <- model$model_info$var_names$longitime
  e1tVar <- "time.prg"
  e2tVar <- "time.trt"
  e0tVar <- "time.cen"
  
  
  # dataset and vectors needed 
  data <- data[data[[longitime]] <= t,] # make sure that there are PSAs
  data.id <- data[!duplicated(data[[idVar]]),]
  data.id <- data.id[data.id[[e2tVar]] >= t, ]
  id_idx <- data.id[[idVar]]
  data <- data[data[[idVar]] %in% id_idx,]
  data.id$t <- t
  data.id$tdt <- t + dt
  n <- nrow(data.id)
  risk_tildet <- rep(NA, n)
  
  # preparation for EPCE
  
  tildT <- pmin(data.id[[e1tVar]], data.id[[e2tVar]], data.id[[e0tVar]], t+dt)
  eid <- rep(NA, n)
  for (i in 1:n) {
    if (data.id[[e1tVar]][i] == min(data.id[[e1tVar]][i], data.id[[e2tVar]][i], data.id[[e0tVar]][i])) {
      eid[i] = 1
    } else if (data.id[[e2tVar]][i] == min(data.id[[e1tVar]][i], data.id[[e2tVar]][i], data.id[[e0tVar]][i])) {
      eid[i] = 2
    } else if (data.id[[e0tVar]][i] == min(data.id[[e1tVar]][i], data.id[[e2tVar]][i], data.id[[e0tVar]][i])) {
      eid[i] = 0
    }
  }
  
  
  osurv_t <- osurv_tildet <- rep(NA, n)
  tilddelta1 <- tilddelta2 <- rep(0,n)
  
  # prediction seed matrix
  set.seed(seed)
  pred_seed <- sample(1:1e6, n)
  
  # Do in parallel 
  ncore <- detectCores() - 1
  plan(multisession, workers = ncore)
  
  out <- lapply(1:n, function(i) {
    future({
      library(rjags)
      ind_data <- data[data[[idVar]] == id_idx[i] & data[[longitime]] <= t,] # & data[[idVar]] here to see if we use only the biomarker until a specific point or the whole 
      visit_vec <- c(t, tildT[i])
      csdypred(model = model, t_start = 1e-6, 
               t_visits = visit_vec, 
               iter = pred_niter, newdata = ind_data, seed = pred_seed[i])
      
    })
  })
  res <- lapply(out, future::value)
  
  for (i in 1:n) {
    pred <- res[[i]]
    est_risk <- pred$risk_visits
    over_surv <- pred$overall.surv.summary
    risk_tildet[i] <- est_risk[est_risk$t.hor == tildT[i], "mean"]
    
    ### indicators for EPCE
    osurv_t[i] <- over_surv[over_surv$t.hor == t, "mean"]
    osurv_tildet[i] <- over_surv[over_surv$t.hor == tildT[i], "mean"]
    if (over_surv[over_surv$t.hor == t, "mean"] == 0) {
      stop("overall survival is 0")
    }
    if (eid[i] == 1 & data.id[[e1tVar]][i] >= t & data.id[[e1tVar]][i] < t + dt) {
      tilddelta1[i] <- 1
    } else {
      tilddelta2[i] <- 1
    }
    
  }
  
  
  data.id$risk_tildet <- risk_tildet
  data.id$tilddelta1 <- tilddelta1
  data.id$tilddelta2 <- tilddelta2
  data.id$osurv_t <- osurv_t
  data.id$osurv_tildet <- osurv_tildet
  
  ### EPCE
  epce_df <- data.frame(id = 1:n,
                        epce = NA)
  epce_df$epce <- -(data.id$tilddelta1 * data.id$risk_tildet + 
      data.id$tilddelta2 * data.id$osurv_tildet / data.id$osurv_t)
  epce <- mean(epce_df$epce)
  
  return(list(data.id = data.id,
              epce_df = epce_df,
              true_epce = epce))
}

# rp <- icjm_auc(model = ICJM1, data = pass[pass$CISNET_ID<=100,], t = 1, dt = 3)
