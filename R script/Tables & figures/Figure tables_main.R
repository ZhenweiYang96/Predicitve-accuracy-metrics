library(tidyverse)
library(ggplot2)
library(cowplot)

# Proportions of each event in simulations --------------------------------
pp <- data.frame(
  id = 1:200,
  prg = NA,
  trt = NA,
  cen = NA
)

for (i in 1:200) {
  load(paste0("Data/Simulation Datasets/traindata_id_", i, ".RData"))
  pp$prg[i] <- sum(train.data.id$status.cmp==1)/nrow(train.data.id)
  pp$trt[i] <- sum(train.data.id$status.cmp==2)/nrow(train.data.id)
  pp$cen[i] <- sum(train.data.id$status.cmp==0)/nrow(train.data.id)
}
round(colMeans(pp),4)

# Comparison between different models -------------------------------------
# AUC and Brier score
ptdf_own <- data.frame(
  dat_id = rep(rep(1:200, each = 3), 4),
  auc = NA,
  bs = NA,
  model = rep(c("Correct", "No baseline", "Linear"), 200 * 4),
  method = rep(factor(c("Naive", "IPCW", "Model-based", "No censoring"),
                      levels = c("Naive", "IPCW", "Model-based", "No censoring")), each = 200 * 3)
)

for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Model-based" &
                 ptdf_own$model == "Correct"] <- res$auc_model
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "IPCW" & 
                 ptdf_own$model == "Correct"] <- res$auc_ipcw
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Model-based" &
                ptdf_own$model == "Correct"] <- res$bs_model
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "IPCW" & 
                ptdf_own$model == "Correct"] <- res$bs_ipcw
  
  load(paste0("Output/Simulation Evaluation/No covariate models/AccMet_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Model-based" &
                 ptdf_own$model == "No baseline"] <- res$auc_model
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "IPCW" & 
                 ptdf_own$model == "No baseline"] <- res$auc_ipcw
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Model-based" &
                ptdf_own$model == "No baseline"] <- res$bs_model
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "IPCW" & 
                ptdf_own$model == "No baseline"] <- res$bs_ipcw
  
  load(paste0("Output/Simulation Evaluation/Linear models/AccMet_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Model-based" &
                 ptdf_own$model == "Linear"] <- res$auc_model
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "IPCW" & 
                 ptdf_own$model == "Linear"] <- res$auc_ipcw
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Model-based" &
                ptdf_own$model == "Linear"] <- res$bs_model
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "IPCW" & 
                ptdf_own$model == "Linear"] <- res$bs_ipcw
  
  # True metrics - no censoring
  load(paste0("Output/Simulation Evaluation/True models/AccMet_true_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "No censoring" &
                 ptdf_own$model == "Correct"] <- true_res$true_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "No censoring" &
                ptdf_own$model == "Correct"] <- true_res$true_bs
  
  load(paste0("Output/Simulation Evaluation/True no cov/AccMet_true_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "No censoring" &
                 ptdf_own$model == "No baseline"] <- true_res$true_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "No censoring" &
                ptdf_own$model == "No baseline"] <- true_res$true_bs
  
  load(paste0("Output/Simulation Evaluation/True linear/AccMet_true_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "No censoring" &
                 ptdf_own$model == "Linear"] <- true_res$true_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "No censoring" &
                ptdf_own$model == "Linear"] <- true_res$true_bs
  
  # Observed metrics - using observed times
  load(paste0("Output/Simulation Evaluation/Observed/Observed models/AccMet_obs_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Naive" &
                 ptdf_own$model == "Correct"] <- obs_res$obs_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Naive" &
                ptdf_own$model == "Correct"] <- obs_res$obs_bs
  
  load(paste0("Output/Simulation Evaluation/Observed/Observed no cov/AccMet_obs_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Naive" &
                 ptdf_own$model == "No baseline"] <- obs_res$obs_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Naive" &
                ptdf_own$model == "No baseline"] <- obs_res$obs_bs
  
  load(paste0("Output/Simulation Evaluation/Observed/Observed linear/AccMet_obs_",
              i, ".RData"))
  ptdf_own$auc[ptdf_own$dat_id == i & 
                 ptdf_own$method == "Naive" &
                 ptdf_own$model == "Linear"] <- obs_res$obs_auc
  ptdf_own$bs[ptdf_own$dat_id == i & 
                ptdf_own$method == "Naive" &
                ptdf_own$model == "Linear"] <- obs_res$obs_bs
}

## AUC
# Plot
auc <- ptdf_own %>% 
  ggplot(aes(y = auc, x = method)) + 
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(aes(yintercept = 0.5), alpha = 0.3) + 
  facet_grid(~ model) + 
  ylab("AUC") +
  xlab("Methods") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 22, face = "bold", family = "LM Roman 10")) # 
auc
ggsave(auc, file="Output/Tables & figures/Main text/AUC_model.png", height = 10, width = 6)

## Brier score
# Plot
bs <- ptdf_own %>% 
  ggplot(aes(y = bs, x = method)) + 
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_hline(aes(yintercept = 0), alpha = 0.3) + 
  facet_grid(~ model) + 
  ylab("Brier Score") +
  xlab("Methods") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 22, face = "bold", family = "LM Roman 10")) # 
bs
ggsave(bs, file="Output/Tables & figures/Main text/BS_model.png", height = 10, width = 6)

## Root mean square error
ptdf_own %>% 
  mutate(auc_true = rep(ptdf_own$auc[ptdf_own$method == "No censoring"], 4),
         bs_true = rep(ptdf_own$bs[ptdf_own$method == "No censoring"], 4)) %>% 
  filter(method != "No censoring") %>% 
  mutate(sqerror_auc = (auc - auc_true)^2,
         sqerror_bs = (bs - bs_true)^2) %>% 
  group_by(model, method) %>% 
  summarise(rmse_auc = sqrt(mean(sqerror_auc)),
            rmse_bs = sqrt(mean(sqerror_bs))) %>% 
  print.data.frame(digits = 2)


## EPCE
ptdf.epce <- data.frame(
  dat_id = rep(rep(1:200, each = 3), 2),
  epce = NA,
  model = rep(c("Correct", "No baseline", "Linear"), 200 * 2),
  method = rep(c("Model-based", "No censoring"), each = 200 * 3)
)

for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "Correct" &
                   ptdf.epce$method == "Model-based"] <- res$epce
  
  load(paste0("Output/Simulation Evaluation/No covariate models/AccMet_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "No baseline" &
                   ptdf.epce$method == "Model-based"] <- res$epce
  
  load(paste0("Output/Simulation Evaluation/Linear models/AccMet_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "Linear" &
                   ptdf.epce$method == "Model-based"] <- res$epce
  
  load(paste0("Output/Simulation Evaluation/True models/EPCE_true_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "Correct" &
                   ptdf.epce$method == "No censoring"] <- true_res$true_epce
  
  load(paste0("Output/Simulation Evaluation/True no cov/EPCE_true_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "No baseline" &
                   ptdf.epce$method == "No censoring"] <- true_res$true_epce
  
  load(paste0("Output/Simulation Evaluation/True linear/EPCE_true_",
              i, ".RData"))
  ptdf.epce$epce[ptdf.epce$dat_id == i & 
                   ptdf.epce$model == "Linear" &
                   ptdf.epce$method == "No censoring"] <- true_res$true_epce
}

epce <- ptdf.epce %>% 
  ggplot(aes(y = epce, x = method)) + 
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_hline(aes(yintercept = 0), alpha = 0.3) + 
  facet_grid(~ model) + 
  ylab("EPCE") +
  xlab("Methods") +
  ylim(c(-0.8, -0.6)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        text = element_text(size = 22, face = "bold", family = "LM Roman 10")) # 

epce
ggsave(epce, file="Output/Tables & figures/Main text/EPCE_model.png", height = 10, width = 6)


# Comparison between schedules --------------------------------------------
# AUC and Brier score
ptdf_sch <- data.frame(
  dat_id = rep(rep(1:200, each = 4), 4), # second is method
  auc = NA,
  bs = NA,
  schedule = rep(factor(c("Random sched \n Unif(0.3, 1)", "Random sched \n Unif(1, 2)", "PASS sched", "Random sched \n Unif(0.3, 4)"), 
                        levels = c("Random sched \n Unif(0.3, 1)", "Random sched \n Unif(1, 2)", "PASS sched", "Random sched \n Unif(0.3, 4)")), 200 * 4),
  method = rep(factor(c("Naive", "IPCW", "Model-based", "No censoring"),
                      levels = c("Naive", "IPCW", "Model-based", "No censoring")), each = 200 * 4)
)

for (i in 1:200) {
  load(paste0("Output/Simulation Evaluation/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Model-based" &
                 ptdf_sch$schedule == "PASS sched"] <- res$auc_model
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "IPCW" & 
                 ptdf_sch$schedule == "PASS sched"] <- res$auc_ipcw
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Model-based" &
                ptdf_sch$schedule == "PASS sched"] <- res$bs_model
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "IPCW" & 
                ptdf_sch$schedule == "PASS sched"] <- res$bs_ipcw
  
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_3_10/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Model-based" &
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- res$auc_model
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "IPCW" & 
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- res$auc_ipcw
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Model-based" &
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- res$bs_model
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "IPCW" & 
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- res$bs_ipcw
  
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_10_20/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Model-based" &
                 ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- res$auc_model
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "IPCW" & 
                 ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- res$auc_ipcw
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Model-based" &
                ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- res$bs_model
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "IPCW" & 
                ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- res$bs_ipcw
  
  load(paste0("Output/Simulation Evaluation/interval_1_4_random_3_40/Correctly specified models/AccMet_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Model-based" &
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- res$auc_model
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "IPCW" & 
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- res$auc_ipcw
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Model-based" &
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- res$bs_model
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "IPCW" & 
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- res$bs_ipcw
  
  # True metrics
  load(paste0("Output/Simulation Evaluation/True models/AccMet_true_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "No censoring"] <- true_res$true_auc
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "No censoring"] <- true_res$true_bs
  
  # Observed metrics - using observed times
  load(paste0("Output/Simulation Evaluation/Observed/Observed models/AccMet_obs_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Naive" &
                 ptdf_sch$schedule == "PASS sched"] <- obs_res$obs_auc
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Naive" &
                ptdf_sch$schedule == "PASS sched"] <- obs_res$obs_bs
  
  load(paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_3_10/AccMet_obs_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Naive" &
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- obs_res$obs_auc
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Naive" &
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 1)"] <- obs_res$obs_bs
  
  load(paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_10_20/AccMet_obs_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Naive" &
                 ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- obs_res$obs_auc
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Naive" &
                ptdf_sch$schedule == "Random sched \n Unif(1, 2)"] <- obs_res$obs_bs
  
  load(paste0("Output/Simulation Evaluation/Observed/Observed interval_1_4_random_3_40/AccMet_obs_",
              i, ".RData"))
  ptdf_sch$auc[ptdf_sch$dat_id == i & 
                 ptdf_sch$method == "Naive" &
                 ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- obs_res$obs_auc
  ptdf_sch$bs[ptdf_sch$dat_id == i & 
                ptdf_sch$method == "Naive" &
                ptdf_sch$schedule == "Random sched \n Unif(0.3, 4)"] <- obs_res$obs_bs
}

## AUC
# Plot
auc <- ptdf_sch %>% 
  ggplot(aes(y = auc, x = method)) + 
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(aes(yintercept = 0.5), alpha = 0.3) + 
  facet_grid(~ schedule) + 
  ylab("AUC") +
  ylim(c(0.2, 0.85))+
  # ggtitle("interval of interest [1, 4)") + 
  xlab("Methods") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 22, face = "bold", family = "LM Roman 10")) # 
auc
ggsave(auc, file="Output/Tables & figures/Main text/AUC_schedule.png", height = 10, width = 10)

## Brier score 
# Plot
bs <- ptdf_sch %>% 
  ggplot(aes(y = bs, x = method)) + 
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_hline(aes(yintercept = 0), alpha = 0.3) + 
  facet_grid(~ schedule) + 
  ylab("Brier Score") +
  # ggtitle("interval of interest [1, 4)") + 
  theme_bw() +
  xlab("Methods") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 22, face = "bold", family = "LM Roman 10"))
bs
ggsave(bs, file="Output/Tables & figures/Main text/BS_schedule.png", height = 10, width = 10)

## Root mean square error
ptdf_sch %>% 
  mutate(auc_true = rep(ptdf_sch$auc[ptdf_sch$method == "No censoring"], 4),
         bs_true = rep(ptdf_sch$bs[ptdf_sch$method == "No censoring"], 4)) %>% 
  filter(method != "No censoring") %>% 
  mutate(sqerror_auc = (auc - auc_true)^2,
         sqerror_bs = (bs - bs_true)^2) %>% 
  group_by(schedule, method) %>% 
  summarise(rmse_auc = sqrt(mean(sqerror_auc)),
            rmse_bs = sqrt(mean(sqerror_bs))) %>% 
  print.data.frame(digits = 2)