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

td_eval <- function(model, data, 
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
  e1tVar <- model$model_info$var_names$survtime[1]
  e2tVar <- model$model_info$var_names$survtime[2]
  Tmax <- max(model$model_info$knots$knot.surv) # longest follow-up in the training set
  
  if (Tmax < t+dt) {
    stop("The interval of interest should be later than the longest follow-up time in the training set")
  }
  
  
  # dataset and vectors needed 
  data <- data[data[[longitime]] <= t,] # make sure that there are PSAs
  data.id <- data[!duplicated(data[[idVar]]),]
  data.id <- data.id[data.id[[e2tVar]] >= t, ]
  id_idx <- data.id[[idVar]]
  data <- data[data[[idVar]] %in% id_idx,]
  data.id$t <- t
  data.id$tdt <- t + dt
  n <- nrow(data.id)
  case_group <- control_group <- 
    w_case_model <- w_control_model <- 
    w_case_ipcw <- w_control_ipcw <- 
    risk_tdt <- rep(NA, n)
  
  # preparation for EPCE
  Qpoints <- matrix(NA, n, qp)
  f.org <- JMbayes:::gaussKronrod(qp)$sk
  w <- JMbayes:::gaussKronrod(qp)$wk
  w <- w[order(f.org)]
  f.org <- f.org[order(f.org)]
  tildT1 <- pmax(data.id[[e1tVar]], t)
  tildT2 <- sapply(1:n, function(i) {
    ifelse(data.id[[eidVar]][i] ==0, t+dt, min(data.id[[e2tVar]][i], t+dt))
  })
  for (i in 1:n) {
    Qpoints[i,] <- (tildT2[i] - tildT1[i])/2 * f.org + (tildT2[i] + tildT1[i])/2
  }
  
  osurv_t <- osurv_tdt <- rep(NA, n)
  tilddelta1 <- tilddelta2 <- rep(0,n)
  risk_qp <- matrix(NA, n, qp)
  
  # IPCW KM preparation
  data.id.km <- data.id
  data.id.km$status.cmp <- ifelse(data.id.km[[eidVar]] == 0, 1, 0)
  km_expression <- formula(paste0("Surv(", e2tVar, ",", eidVar, ")~1"))
  ipcw_km <- survfit(km_expression, data = data.id.km)
  
  # prediction seed matrix
  set.seed(seed)
  pred_seed <- sample(1:1e6, n)
  
  ### Note: we have to make sure that there are absolute cases and controls!
  if (sum(data.id[[e1tVar]] >= t & data.id[[e2tVar]] <= t+dt & data.id[[eidVar]] ==1) == 0) {
    stop(paste0("The IPCW approach needs at least one absolute case. The current interval of interests (",
                t, ", ", t+dt, "] does not contain absolute cases."))
  }
  if (sum(data.id[[e1tVar]] >= t+dt) == 0) {
    stop(paste0("The IPCW approach needs at least one absolute control who survived the current interval of interests (",
                t, ", ", t+dt, "]."))
  }
  
  # Do in parallel 
  ncore <- detectCores() - 1
  plan(multisession, workers = ncore)
  
  out <- lapply(1:n, function(i) {
    future({
      library(rjags)
      ind_data <- data[data[[idVar]] == id_idx[i] & data[[longitime]] <= t,] # & data[[idVar]] here to see if we use only the biomarker until a specific point or the whole 
      if (!(data.id[[e1tVar]][i] <= t+dt  & data.id[[e2tVar]][i] >= t)) {
        visit_vec <- c(t, t+dt, Tmax)
      } else {
        visit_vec <- c(t, t+dt, 
                       ifelse(data.id[[e1tVar]][i] == 0, 1e-5, data.id[[e1tVar]][i]), 
                       data.id[[e2tVar]][i], 
                       Qpoints[i,], Tmax)
      }
      csdypred(model = model, t_start = 1e-6, 
               t_visits = visit_vec, 
               iter = pred_niter, newdata = ind_data, seed = pred_seed[i], trt_part = T)
      
    })
  })
  res <- lapply(out, future::value)
  
  for (i in 1:n) {
    
    #### Weights for Time-dependent AUC and Brier score
    pred <- res[[i]]
    est_risk <- pred$risk_visits
    over_surv <- pred$overall.surv.summary
    risk_tdt[i] <- est_risk[est_risk$t.hor == t+dt, "mean"]
    # special case for left time point == 0, let it to be 1e-5
    time.left <- ifelse(data.id[[e1tVar]][i] == 0, 1e-5, data.id[[e1tVar]][i])
    if (data.id[[e1tVar]][i] <= t+dt) {
      LFTCRT_RISK <- est_risk[est_risk$t.hor == time.left, "mean"]
      LFTCRT_SURV <- over_surv[over_surv$t.hor == time.left, "mean"]
    }
    
    if (data.id[[e2tVar]][i] <= t) {
      case_group[i] <- 0
      w_case_model[i] <- 0
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 0
      w_control_model[i] <- 0
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] <= t & data.id[[e2tVar]][i] >= t & data.id[[e2tVar]][i] <= t+dt & data.id[[eidVar]][i] == 0) {
      case_group[i] <- 1
      w_case_model[i] <- (est_risk[est_risk$t.hor == t+dt, "mean"] - est_risk[est_risk$t.hor == t, "mean"])/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 1
      w_control_model[i] <- est_risk[est_risk$t.hor == Inf, "mean"] - 
        (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] <= t & data.id[[e2tVar]][i] >= t & data.id[[e2tVar]][i] <= t+dt & data.id[[eidVar]][i] != 0) {
      case_group[i] <- 2
      w_case_model[i] <- (est_risk[est_risk$t.hor == data.id[[e2tVar]][i], "mean"] - est_risk[est_risk$t.hor == t, "mean"])/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 0
      w_control_model[i] <- 0
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] >= t & data.id[[e1tVar]][i] <= t+dt & data.id[[e2tVar]][i] >= t & data.id[[e2tVar]][i] <= t+dt & data.id[[eidVar]][i] == 0) {
      case_group[i] <- 3
      w_case_model[i] <- (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 1
      w_control_model[i] <- est_risk[est_risk$t.hor == Inf, "mean"] - 
        (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] >= t & data.id[[e1tVar]][i] <= t+dt & data.id[[e2tVar]][i] >= t & data.id[[e2tVar]][i] <= t+dt & data.id[[eidVar]][i] == 1) {
      case_group[i] <- 4
      w_case_model[i] <- 1
      est_risk_ipcw <- summary(ipcw_km, times = c(t, data.id[[e2tVar]][i]))
      w_case_ipcw[i] <- 1/(est_risk_ipcw$surv[2]/est_risk_ipcw$surv[1])
      
      control_group[i] <- 0
      w_control_model[i] <- 0
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] >= t & data.id[[e1tVar]][i] <= t+dt & data.id[[e2tVar]][i] >= t & data.id[[e2tVar]][i] <= t+dt & data.id[[eidVar]][i] == 2) {
      case_group[i] <- 5
      w_case_model[i] <- (est_risk[est_risk$t.hor == data.id[[e2tVar]][i], "mean"] - LFTCRT_RISK)/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 0
      w_control_model[i] <- 0
      w_control_ipcw[i] <- 0
    } else if (data.id[[e1tVar]][i] >= t & data.id[[e1tVar]][i] <= t+dt & data.id[[e2tVar]][i] >= t+dt) {
      case_group[i] <- 3
      w_case_model[i] <- (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      if (data.id[[eidVar]][i] == 0) {
        control_group[i] <- 1
        w_control_model[i] <- est_risk[est_risk$t.hor == Inf, "mean"] - 
          (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
        w_control_ipcw[i] <- 0
      } else {
        control_group[i] <- 2
        w_control_model[i] <- (est_risk[est_risk$t.hor == data.id[[e2tVar]][i], "mean"] - est_risk[est_risk$t.hor == t+dt, "mean"])/LFTCRT_SURV
        w_control_ipcw[i] <- 0
      }
    } else if (data.id[[e1tVar]][i] <= t & data.id[[e2tVar]][i] >= t+dt) {
      case_group[i] <- 1
      w_case_model[i] <- (est_risk[est_risk$t.hor == t+dt, "mean"] - est_risk[est_risk$t.hor == t, "mean"])/LFTCRT_SURV
      w_case_ipcw[i] <- 0
      
      if (data.id[[eidVar]][i] == 0) {
        control_group[i] <- 1
        w_control_model[i] <- est_risk[est_risk$t.hor == Inf, "mean"] - 
          (est_risk[est_risk$t.hor == t+dt, "mean"] - LFTCRT_RISK)/LFTCRT_SURV
        w_control_ipcw[i] <- 0
      } else {
        control_group[i] <- 2
        w_control_model[i] <- (est_risk[est_risk$t.hor == data.id[[e2tVar]][i], "mean"] - est_risk[est_risk$t.hor == t+dt, "mean"])/LFTCRT_SURV
        w_control_ipcw[i] <- 0
      }
    } else if (data.id[[e2tVar]][i] >= t+dt) {
      case_group[i] <- 0
      w_case_model[i] <- 0
      w_case_ipcw[i] <- 0
      
      control_group[i] <- 3
      w_control_model[i] <- 0
      est_risk_ipcw <- summary(ipcw_km, times = c(t, t+dt))
      w_control_ipcw[i] <- 1/(est_risk_ipcw$surv[2]/est_risk_ipcw$surv[1])
    }
    
    ### indicators for EPCE
    osurv_t[i] <- over_surv[over_surv$t.hor == t, "mean"]
    if (over_surv[over_surv$t.hor == t, "mean"] == 0) {
      stop("overall survival is 0")
    }
    osurv_tdt[i] <- over_surv[over_surv$t.hor == t+dt, "mean"]
    
    colnames(risk_qp) <- sapply(1:qp, function(x) {paste0("qprisk",x)})
    if (data.id[[e1tVar]][i] <= t+dt  & data.id[[e2tVar]][i] >= t) {
      tilddelta1[i] <- 1
      risk_qp[i,] <- est_risk[est_risk$t.hor %in% Qpoints[i,], "mean"]
    } else {
      risk_qp[i,] <- 0
    }
    if (data.id[[e1tVar]][i] > t+dt | data.id[[eidVar]][i] == 0) {
      tilddelta2[i] <- 1
    }
    
  }
  
  data.id$case_group <- case_group
  data.id$control_group <- control_group
  data.id$w_case_model <- w_case_model
  data.id$w_case_ipcw <- w_case_ipcw
  data.id$w_control_model <- w_control_model
  data.id$w_control_ipcw <- w_control_ipcw
  data.id$risk_tdt <- risk_tdt
  data.id$tilddelta1 <- tilddelta1
  data.id$tilddelta2 <- tilddelta2
  data.id$osurv_t <- osurv_t
  data.id$osurv_tdt <- osurv_tdt
  data.id <- cbind(data.id, risk_qp)
  
  
  ### Time-dependent AUC
  auc_df <- data.frame(c = seq(0, 1, length.out=c.step),
                       sens_model = NA,
                       speinvs_model = NA,
                       sens_ipcw = NA,
                       speinvs_ipcw = NA)
  
  for (i in 1:nrow(auc_df)) {
    auc_df$sens_model[i] <- sum((data.id$risk_tdt >= auc_df$c[i]) * data.id$w_case_model)/sum(data.id$w_case_model)
    auc_df$speinvs_model[i] <- 1 - sum((data.id$risk_tdt < auc_df$c[i]) * data.id$w_control_model)/sum(data.id$w_control_model)
    auc_df$sens_ipcw[i] <- sum((data.id$risk_tdt >= auc_df$c[i]) * data.id$w_case_ipcw)/sum(data.id$w_case_ipcw)
    auc_df$speinvs_ipcw[i] <- 1 - sum((data.id$risk_tdt < auc_df$c[i]) * data.id$w_control_ipcw)/sum(data.id$w_control_ipcw)
  }
  
  auc_model_df <- auc_df[!duplicated(auc_df$speinvs_model),]
  auc_model <- integrate(approxfun(auc_model_df$speinvs_model, auc_model_df$sens_model), lower = 0, upper = 1)
  auc_ipcw_df <- auc_df[!duplicated(auc_df$speinvs_ipcw),]
  auc_ipcw <- integrate(approxfun(auc_ipcw_df$speinvs_ipcw, auc_ipcw_df$sens_ipcw), lower = 0, upper = 1)
  
  ### Time-dependent Brier score
  bs_df <- data.frame(id = 1:n,
                      bs_model = NA,
                      bs_ipcw = NA)
  bs_df$bs_model <- (1 - data.id$risk_tdt)^2 * data.id$w_case_model + (0 - data.id$risk_tdt)^2 * data.id$w_control_model
  bs_df$bs_ipcw <- (1 - data.id$risk_tdt)^2 * data.id$w_case_ipcw + (0 - data.id$risk_tdt)^2 * data.id$w_control_ipcw
  bs_model <- mean(bs_df$bs_model)
  bs_ipcw <- mean(bs_df$bs_ipcw)
  
  ### EPCE
  epce_df <- data.frame(id = 1:n,
                        epce = NA)
  epce_df$epce <- - (data.id$tilddelta1 * (risk_qp %*% w) + 
      data.id$tilddelta2 * data.id$osurv_tdt / data.id$osurv_t)
  epce <- mean(epce_df$epce)
  
  return(list(data.id = data.id,
              auc_df = auc_df,
              auc_model = auc_model$value,
              auc_ipcw = auc_ipcw$value,
              bs_df = bs_df,
              bs_model = bs_model,
              bs_ipcw = bs_ipcw,
              epce_df = epce_df,
              epce = epce))
}

# rp <- icjm_auc(model = ICJM1, data = pass[pass$CISNET_ID<=100,], t = 1, dt = 3)
