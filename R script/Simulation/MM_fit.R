# Package -----------------------------------------------------------------
library(splines)
library(GLMMadaptive)
library(future)
library(parallel)

# Function ----------------------------------------------------------------
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
out <- lapply(c(1:200), function(datanum) {
  future({
    #Real simulation/True sens 75/
    load(paste0("Data/Simulation Datasets/traindata_",
                datanum, ".RData"))
    
    # train.data <- train.data[train.data$CISNET_ID %in% 1:100,]
    # train.biopsy <- train.biopsy[train.biopsy$CISNET_ID %in% 1:100,]
    
    # start modelling
    mm <- mm_fit(data = train.data)
    save(mm,
         file = paste0("Output/Simulation Models/Fitted mixed model/MM_",
                       datanum, ".RData"))
    #
    plan(multisession, workers = ncore)
    
  })
})
res <- lapply(out, future::value)


# Extract inits from the mixed models -------------------------------------
for (datanum in 1:200) {
  load(paste0("Output/Simulation Models/Fitted mixed model/MM_",
              datanum, ".RData"))
  mm_inits <- list(betaL = fixef(mm),
                   b = ranef(mm),
                   inv_D = solve(mm$D),
                   tau = 1/exp(mm$phis)^2)
  save(mm_inits,
       file = paste0("Output/Simulation Models/Fitted mixed model/MM_inits_",
                     datanum, ".RData"))
}

