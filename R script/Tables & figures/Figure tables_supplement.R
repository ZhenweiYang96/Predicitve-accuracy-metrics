library(tidyverse)
library(ggplot2)
library(cowplot)

load("Output/Simulation Models/Correctly specified models/Joint model_1.RData")
knot.longi <- icjm$model_info$knots$knot.longi
load("Output/Simulation Storage/corspe/pred_1.RData")
pred1 <- res
load("Output/Simulation Storage/linear/pred_1.RData")
pred2 <- res
load("Data/Simulation Datasets/traindata_2.RData")

set.seed(2025)
train.data <- train.data[train.data$TimeSince_Dx <= 1 & train.data$time.cmp2 >= 1,]
train.data.id <- train.data[!duplicated(train.data$CISNET_ID),]
samples <- sample(1:nrow(train.data.id), 8)
samples_id <- train.data.id$CISNET_ID[samples]
lent <- 2
sp <- ns(seq(0, lent, length.out=50), knots = knot.longi[2:3], B = knot.longi[c(1,4)])


psa_true <- data.frame(
  id = rep(rep(samples_id, each=50), 2),
  time = rep(rep(seq(0, lent, length.out=50), length(samples_id)), 2),
  time1 = rep(rep(sp[,1], length(samples_id)), 2),
  time2 = rep(rep(sp[,2], length(samples_id)), 2),
  time3 = rep(rep(sp[,3], length(samples_id)), 2),
  age = rep(rep(train.data.id$DxAge[samples], each=50), 2),
  Model = rep(c("Correct", "Linear"), each = 50 * length(samples_id)),
  psa = NA,
  psa_lower = NA,
  psa_upper = NA
)

psa_observed <- train.data[train.data$CISNET_ID %in% samples_id,] %>% 
  mutate(id = CISNET_ID)


for (i in 1:length(samples)){
  para <- pred1[[samples[i]]]$parameters
  beta <- para$beta
  b <- para$b
  age <- train.data.id$DxAge[train.data.id$CISNET_ID == samples_id[i]]
  psa_true$psa[psa_true$id == samples_id[i] & psa_true$Model == "Correct"] <- sapply(seq(0, lent, length.out=50), function(n) {
    mean(sapply(1:nrow(beta), function(l) {
      c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age) %*% beta[l,] +
        c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)])) %*% b[l,]
    }))
  })
  psa_true$psa_lower[psa_true$id == samples_id[i] & psa_true$Model == "Correct"] <- sapply(seq(0, 3, length.out=50), function(n) {
    quantile(sapply(1:nrow(beta), function(l) {
      c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age) %*% beta[l,] +
        c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)])) %*% b[l,]
    }), probs = 0.025)
  })
  psa_true$psa_upper[psa_true$id == samples_id[i] & psa_true$Model == "Correct"] <- sapply(seq(0, 3, length.out=50), function(n) {
    quantile(sapply(1:nrow(beta), function(l) {
      c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age) %*% beta[l,] +
        c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)])) %*% b[l,]
    }), probs = 0.975)
  })
  
  # Linear model
  para <- pred2[[samples[i]]]$parameters
  beta <- para$beta
  b <- para$b
  age <- train.data.id$DxAge[train.data.id$CISNET_ID == samples_id[i]]
  psa_true$psa[psa_true$id == samples_id[i] & psa_true$Model == "Linear"] <- sapply(seq(0, lent, length.out=50), function(n) {
    mean(sapply(1:nrow(beta), function(l) {
      c(1, n, age) %*% beta[l,] +
        c(1, n) %*% b[l,]
    }))
  })
  psa_true$psa_lower[psa_true$id == samples_id[i] & psa_true$Model == "Linear"] <- sapply(seq(0, lent, length.out=50), function(n) {
    quantile(sapply(1:nrow(beta), function(l) {
      c(1, n, age) %*% beta[l,] +
        c(1, n) %*% b[l,]
    }), probs = 0.025)
  })
  psa_true$psa_upper[psa_true$id == samples_id[i] & psa_true$Model == "Linear"] <- sapply(seq(0, lent, length.out=50), function(n) {
    quantile(sapply(1:nrow(beta), function(l) {
      c(1, n, age) %*% beta[l,] +
        c(1, n) %*% b[l,]
    }), probs = 0.975)
  })
  
}


# plot
psa <- ggplot(psa_true) +
  geom_line(data = psa_true %>% filter(time <=1), aes(x = time, y = psa, color = Model)) + 
  geom_line(data = psa_true %>% filter(time >1), aes(x = time, y = psa, color = Model), linetype="dashed") + 
  geom_ribbon(aes(x = time, ymin = psa_lower, ymax = psa_upper, fill = Model), alpha = 0.3) +
  geom_point(data = psa_observed, aes(x = TimeSince_Dx, y = PSAValue)) + 
  facet_wrap(~ id, nrow = 2) + 
  ylab(expression(log[2](PSA + 1))) +
  theme_bw() +
  scale_x_continuous(
    name = "Time (years)",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2),     # Set custom breaks
    labels = c("0", "0.5", "1", "1.5", "2")  # Custom labels
  ) + 
  theme(legend.position = "top",
        text = element_text(size = 22, face = "bold", family = "LM Roman 10"))
psa
ggsave(psa, file="Output/Tables & figures/Supplementary/PSAtrajecotry.png", height = 8, width = 10)
