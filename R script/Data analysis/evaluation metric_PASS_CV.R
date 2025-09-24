# Split the data ----------------------------------------------------------
library(tidyverse)
set.seed(100)
load("Data/pass_id.RData")
load("Data/pass.RData")
length(pass.id$CISNET_ID)
od.id <- sample(833,833)
pass.id$newod <- od.id
pass.id <- pass.id %>% 
  mutate(group_id = case_when(newod %in% 1:166 ~ 1,
                              newod %in% 167:332 ~ 2,
                              newod %in% 333:498 ~ 3,
                              newod %in% 499:664 ~ 4,
                              newod %in% 665:833 ~ 5))
gid1 <- pass.id$CISNET_ID[pass.id$group_id == 1]
gid2 <- pass.id$CISNET_ID[pass.id$group_id == 2]
gid3 <- pass.id$CISNET_ID[pass.id$group_id == 3]
gid4 <- pass.id$CISNET_ID[pass.id$group_id == 4]
gid5 <- pass.id$CISNET_ID[pass.id$group_id == 5]

# Group 1
train.pass.id <- pass.id[!(pass.id$CISNET_ID %in% gid1),]
train.pass <- pass[!(pass$CISNET_ID %in% gid1),]
save(train.pass.id, file="Data/CV/train_pass_id_1.RData")
save(train.pass, file="Data/CV/train_pass_1.RData")

test.pass.id <- pass.id[pass.id$CISNET_ID %in% gid1,]
test.pass <- pass[pass$CISNET_ID %in% gid1,]
save(test.pass.id, file="Data/CV/test_pass_id_1.RData")
save(test.pass, file="Data/CV/test_pass_1.RData")

# Group 2
train.pass.id <- pass.id[!(pass.id$CISNET_ID %in% gid2),]
train.pass <- pass[!(pass$CISNET_ID %in% gid2),]
save(train.pass.id, file="Data/CV/train_pass_id_2.RData")
save(train.pass, file="Data/CV/train_pass_2.RData")

test.pass.id <- pass.id[pass.id$CISNET_ID %in% gid2,]
test.pass <- pass[pass$CISNET_ID %in% gid2,]
save(test.pass.id, file="Data/CV/test_pass_id_2.RData")
save(test.pass, file="Data/CV/test_pass_2.RData")

# Group 3
train.pass.id <- pass.id[!(pass.id$CISNET_ID %in% gid3),]
train.pass <- pass[!(pass$CISNET_ID %in% gid3),]
save(train.pass.id, file="Data/CV/train_pass_id_3.RData")
save(train.pass, file="Data/CV/train_pass_3.RData")

test.pass.id <- pass.id[pass.id$CISNET_ID %in% gid3,]
test.pass <- pass[pass$CISNET_ID %in% gid3,]
save(test.pass.id, file="Data/CV/test_pass_id_3.RData")
save(test.pass, file="Data/CV/test_pass_3.RData")

# Group 4
train.pass.id <- pass.id[!(pass.id$CISNET_ID %in% gid4),]
train.pass <- pass[!(pass$CISNET_ID %in% gid4),]
save(train.pass.id, file="Data/CV/train_pass_id_4.RData")
save(train.pass, file="Data/CV/train_pass_4.RData")

test.pass.id <- pass.id[pass.id$CISNET_ID %in% gid4,]
test.pass <- pass[pass$CISNET_ID %in% gid4,]
save(test.pass.id, file="Data/CV/test_pass_id_4.RData")
save(test.pass, file="Data/CV/test_pass_4.RData")

# Group 5
train.pass.id <- pass.id[!(pass.id$CISNET_ID %in% gid5),]
train.pass <- pass[!(pass$CISNET_ID %in% gid5),]
save(train.pass.id, file="Data/CV/train_pass_id_5.RData")
save(train.pass, file="Data/CV/train_pass_5.RData")

test.pass.id <- pass.id[pass.id$CISNET_ID %in% gid5,]
test.pass <- pass[pass$CISNET_ID %in% gid5,]
save(test.pass.id, file="Data/CV/test_pass_id_5.RData")
save(test.pass, file="Data/CV/test_pass_5.RData")


# Fit the models ----------------------------------------------------------

### MM fit
library(splines)
library(GLMMadaptive)
library(future)
library(parallel)

mm_fit <- function(df = 3, data) {
  
  # Start with the mixed model ----------------------------------------------
  bound <- range(data$TimeSince_Dx, data$time.cmp1, data$time.cmp2)
  data[,(ncol(data)+1):(ncol(data)+3)] <- ns(data$TimeSince_Dx, 3, B = bound)
  colnames(data)[(ncol(data)-2):ncol(data)] <- c("TimeSince_Dx.1",
                                                 "TimeSince_Dx.2",
                                                 "TimeSince_Dx.3")
  mm <- mixed_model(fixed = PSAValue ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 + DxAge,
                    data = data, 
                    random = ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 | CISNET_ID,
                    family = students.t(df = df), n_phis = 1,
                    initial_values = list("betas" = gaussian()))
}


# Fit the mixed models ----------------------------------------------------
ncore <- detectCores() - 1
plan(multisession, workers = ncore)
out <- lapply(c(1:5), function(datanum) {
  future({
    #Real simulation/True sens 75/
    load(paste0("Data/CV/train_pass_",
                datanum, ".RData"))
    
    # start modelling
    mm <- mm_fit(data = train.pass)
    save(mm,
         file = paste0("Data/CV/Fitted mixed model/MM_",
                       datanum, ".RData"))

    
  })
})
res <- lapply(out, future::value)

for (datanum in 1:5) {
  load(paste0("Data/CV/Fitted mixed model/MM_",
              datanum, ".RData"))
  mm_inits <- list(betaL = fixef(mm),
                   b = ranef(mm),
                   inv_D = solve(mm$D),
                   tau = 1/exp(mm$phis)^2)
  save(mm_inits,
       file = paste0("Data/CV/Fitted mixed model/MM_inits_",
                     datanum, ".RData"))
}


# Fit the model -----------------------------------------------------------
library(future)
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
seed = 2008:2012

# correctly-specified models -------------------------------------
source("R script/Function/ICJM_correctlyspecify.R")
ncore_div3 <- floor(ncore/3)
plan(
  list(
    tweak(multisession, workers = ncore_div3),
    tweak(multisession, workers = 3)
  )
)

t1 <- Sys.time()
out <- lapply(1:ncore_div3, function(datanum) {
  future({
    
    # load data file
    load(paste0("Data/CV/train_pass_",
                datanum, ".RData"))
    # load inits
    load(paste0("Data/CV/Fitted mixed model/MM_inits_",
                datanum, ".RData"))
    
    # start modelling
    icjm <- icjm_cor(data = train.pass, n.adapt = 10000, 
                     n.burnin = 10000, n.iter = 10000, 
                     GLMM_fit = mm_inits,
                     seed = seed[datanum])
    save(icjm,
         file = paste0("Data/CV/ICJM_CV_",
                       datanum, ".RData"))
    # 
    plan(
      list(
        tweak(multisession, workers = ncore_div3),
        tweak(multisession, workers = 3)
      )
    )
    return(paste0(datanum,"done!"))
  })
})
res <- lapply(out, future::value)
t2 <- Sys.time()

t2 - t1


# Calculating the Metrics -------------------------------------------------
source("R script/Function/time dependent evaluation metrics.R")
seed = 2008:2012
for (i in 1:5) {
  load(paste0("Data/CV/ICJM_CV_", i, ".RData"))
  load(paste0("Data/CV/test_pass_",i , ".RData"))
  res <- td_eval(model = icjm, data = test.pass, type = "corspe", 
                 t = 1, dt = 3, seed = seed[i])
  save(res, file=paste0("Output/Data analysis/AccMet_CV_", i, ".RData"))
}

# Output ------------------------------------------------------------------
eva_pass <- data.frame(
  cvround = 1:5,
  auc_model = NA,
  auc_ipcw = NA,
  bs_model = NA,
  bs_ipcw = NA,
  epce = NA
)

for (i in 1:5) {
  load(paste0("Output/Data analysis/AccMet_CV_", i, ".RData"))
  eva_pass$auc_model[i] <- res$auc_model
  eva_pass$auc_ipcw[i] <- res$auc_ipcw
  eva_pass$bs_model[i] <- res$bs_model
  eva_pass$bs_ipcw[i] <- res$bs_ipcw
  eva_pass$epce[i] <- res$epce
}
colMeans(eva_pass)

