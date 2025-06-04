# # package
library(splines)
library(mvtnorm)
library(MASS)
library(rjags)
source("Function/function.R")

# risk prediction function  -----------------------------------------------

csdypred <- function(model = NULL,
                     t_start = 0,
                     t_visits = NULL,
                     n.step = 50L,
                     iter = 250L, n.adapt = 150L,
                     seed = 100L,
                     idVar = "CISNET_ID",
                     varname = list(longitimeVar = NULL,
                                    #survtimeVar = NULL,
                                    longi.bsVar = NULL,
                                    surv.bsVar = NULL),
                     newdata,
                     sigma.inits = 0.27,
                     sigMat.inits = NULL,
                     qp = 15,
                     trt_part = F) {
  if (is.null(newdata)) {
    stop("Test set is missing! \n")
  }
  
  if (is.null(t_visits)) {
    stop("The clinical visits time points are missing\n")
  }
  
  if (length(unique(newdata[[idVar]])) > 1) {
    stop("Dynamic predictions can be generated for only one patient at a time.\n")
  }
  
  if (is.null(model)) {
    stop("model should be provided")
  } else {
    object <- model
  }
  
  if (min(t_visits) < t_start) {
    stop("start point should be not be later than visit times!")
  } 
  
  Survtimes <- sort(c(t_visits,
                      seq(t_start, max(t_visits), length.out = n.step)))
  Survtimes <- Survtimes[!duplicated(Survtimes)] # remove duplicated time points
  
  u <- max(Survtimes)
  t <- min(Survtimes)
  
  # extract the posterior distribution
  
  jags_res <- extract_pd(object$mcmc)
  pool.bs_gammas <- jags_res$pool.bs_gammas
  pool.gammas <- jags_res$pool.gammas
  pool.alphas <- jags_res$pool.alphas
  pool.beta <- jags_res$pool.beta
  pool.D <- jags_res$pool.D
  pool.tau <- jags_res$pool.tau
  
  # extract model info from JAGS
  # idVar <- ifelse(is.null(varname$idVar),
  #                 object$model_info$var_names$id,
  #                 varname$idVar)
  
  test <- newdata
  test.id <- test[!duplicated(test[[idVar]]),]
  df.ns.longi <- sqrt(ncol(pool.D)) - 1
  knot.longi <- object$model_info$knots$knot.longi
  knot.surv <- object$model_info$knots$knot.surv
  
  
  # survtimeVar <- ifelse(is.null(varname$survtimeVar),
  #                       object$model_info$var_names$survtime,
  #                       varname$survtimeVar)
  longitimeVar <- ifelse(is.null(varname$longitimeVar),
                         object$model_info$var_names$longitime,
                         varname$longitimeVar)
  longi.bsVar <- ifelse(is.null(varname$longi.bsVar),
                        object$model_info$var_names$longi.bsVar,
                        varname$longi.bsVar)
  surv.bsVar <- ifelse(is.null(varname$surv.bsVar),
                       object$model_info$var_names$surv.bsVar,
                       varname$surv.bsVar)
  
  # preparation for sample from posterior distribution
  step <- Survtimes
  total <- iter
  
  beta <- matrix(NA, total, ncol(pool.beta))
  b <- matrix(NA, total, df.ns.longi + 1)
  tau <- matrix(NA, total, 1)
  gam_bh <- array(NA, dim = c(total, ncol(pool.bs_gammas)/2, 2))
  gamma <- array(NA, dim =c(total, ncol(pool.gammas)/2, 2))
  Sigma <- array(NA, dim = c(total, df.ns.longi + 1, df.ns.longi + 1))
  alpha <- array(NA, dim = c(total, ncol(pool.alphas)/2, 2)) # col = risk
  
  # risk output data frame
  risk <- data.frame(t.hor = step,
                     mean = rep(NA, length(step)),
                     upper = rep(NA, length(step)),
                     lower = rep(NA, length(step)))
  overall.surv.summary <- data.frame(t.hor = step,
                                     mean = rep(NA, length(step)),
                                     upper = rep(NA, length(step)),
                                     lower = rep(NA, length(step)))
  # fitted longitudinal model
  longi <- data.frame(t.hor = seq(0, max(test[[longitimeVar]]), length.out = 50),
                      PSA.pred = rep(NA,50))
  
  f.org <- JMbayes:::gaussKronrod(qp)$sk
  w <- JMbayes:::gaussKronrod(qp)$wk
  
  # risk store set
  risk.store <- matrix(NA, iter, length(step))
  
  # sample
  set.seed(seed)
  idx <- sample(1:nrow(pool.beta), total, replace = T)
  beta <-  pool.beta[idx,]
  
  tau[,1] <- pool.tau[idx]
  
  for (i in 1:ncol(pool.bs_gammas)/2) {
    for (j in 1:2) {
      gam_bh[,i,j] <- pool.bs_gammas[idx, i + (j-1) * ncol(pool.bs_gammas)/2]
    }
  }
  
  for (i in 1:ncol(pool.gammas)/2) {
    for (j in 1:2) {
      gamma[,i,j] <- pool.gammas[idx, i + (j-1) * ncol(pool.gammas)/2]
    }
  }
  
  for (i in 1:length(idx)) {
    Sigma[i,,] <- pool.D[idx[i],]
  }
  
  for (i in 1:length(idx)) {
    alpha[i,,] <- pool.alphas[idx[i],]
  }
  
  #### sample random effects
  n_test <- nrow(test)
  XL <- cbind(1,
              ns(test[[longitimeVar]], knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]),
              test[[longi.bsVar]])
  ZL <- cbind(1,
              ns(test[[longitimeVar]], knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
  XL.b.now <- matrix(1, qp, ncol(pool.beta))
  ZL.b.now <- matrix(1, qp, df.ns.longi+1)
  XL.b.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  XL.b.now[,2:(df.ns.longi+1)] <- ZL.b.now[,2:(df.ns.longi+1)] <- ns(f.org * t/2 + t/2,
                                                                     knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  
  XL.b.back <- matrix(1, qp, ncol(pool.beta))
  ZL.b.back <- matrix(1, qp, df.ns.longi+1)
  XL.b.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  XL.b.back[,2:(df.ns.longi+1)] <- ZL.b.back[,2:(df.ns.longi+1)] <- ns(f.org * t/2 + t/2 - 1,
                                                                       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  #outer(ns(test.id$time, knots = knot.longi[2:3], B = knot.longi[c(1,4)])/2, f.org + 1)
  
  X <- test.id[[surv.bsVar]]
  T.b <- splineDesign(knot.surv, outer(t/2, f.org + 1), ord = 4L)
  
  
  loglike <- function(b) {
    n_test/2 * log(tau[i,] / (3*pi))  - 2 * sum(sapply(1:n_test, function(n) {
      log(1 + tau[i,]/3 * (test$PSAValue[n] - (XL[n,,drop = F] %*% beta[i,] + ZL[n,, drop = F] %*% b))^2)
    }))  -
      # n_test/2 * log(tau[i,]) - 0.5 * tau[i,] * sum((test$PSAValue - (XL %*% beta[i,] + ZL %*% b))^2) -
      t/2 * sum(sapply(1:2, function(k) {
        exp(gamma[i,,k] * X) * w %*% sapply(1:qp, function(l) {
          exp(T.b[l,] %*% gam_bh[i,,k] + alpha[i,1,k] * (XL.b.now[l,] %*% beta[i,] + ZL.b.now[l,] %*% b) +
                alpha[i,2,k] * (XL.b.now[l,] %*% beta[i,] + ZL.b.now[l,] %*% b -
                                  XL.b.back[l,] %*% beta[i,] - ZL.b.back[l,]%*% b))
        })
      })) -
      0.5 * log(det(Sigma[i,,])) - 0.5 * matrix(b,1,4) %*% solve(Sigma[i,,]) %*% matrix(b,4,1)
  }
  
  # RWMH method -------------------------------------------------------------
  update.sigma<- function(sigma2, acc, p=p, j) {
    c=((1-1/d)*sqrt(2*pi)*exp(a^2/2)/(2*a) + 1/(d*p*(1-p)))
    Theta=log(sqrt(sigma2))
    Theta=Theta+c*(acc-p)/max(200, i/d)
    return(exp(Theta))
  }
  
  
  update.cov<-function(sigMat, j, thetaM, theta){
    #function to recursively update covariance matrix, as part of adaptive MCMC sampling, updating the covariance matrix
    epsilon=1/j
    thetaM2=((thetaM*j)+theta)/(j+1)
    sigMat=(j-1)/j*sigMat + thetaM%*%t(thetaM)-(j+1)/j*thetaM2%*%t(thetaM2)+1/j*theta%*%t(theta) + epsilon*diag(d)
    return(list(sigMat=sigMat, thetaM=thetaM2))
  }
  
  
  #############################################
  # Begin
  #initialise the RM section
  #############################################
  #the optimal acceptance probability for the multivariate case
  pstar=0.234
  a=-qnorm(pstar/2)
  n0=round(5/(pstar*(1-pstar)))
  #iMax, is the max number of iterations before the last restart
  iMax=20
  Numbig=0
  Numsmall=0
  
  #initialise the MCMC program
  #niter=iter + n.adapt      # niter is the number of iterations in the Markov chain.
  # Change it to the value you want.
  d=df.ns.longi+1 #dimension of parameters to be updated
  output<-rep(0, d) #output parameter values
  output.mat<-output #records all MCMC output for the parameters
  # initial values
  sigma <- sigma.inits
  if (is.null(sigMat.inits)) {sigMat <- diag(d)} else {sigMat <- sigMat.inits}
  meanacc.record <- list()
  # end ---------------------------------------------------------------------
  
  #### End adaption and start MH
  n.adapt <- c(n.adapt, rep(30, iter - 1))
  for (i in 1:iter) {
    
    #############################################
    # Begin
    #initialise the RM section
    #############################################
    sigma= sigma       #an arbitrary starting value for sigma, equivalent to theta=ln(sigma)=0
    sigma2=sigma^2 #variance
    sigma.start<- sigma
    sigma.vec<- sigma
    sigMat <- sigMat
    acc.vec<-rep(NA, n.adapt[i]) #records the acceptance probabilities from each MCMC iteration
    num.restart=0
    j=1
    #############################################
    # End
    #initialise the RM section
    #############################################
    
    for (m in 2:n.adapt[i]) {
      
      #propose a new value of theta
      output.prop<- as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
      pi.old.log <- loglike(output)
      pi.new.log <- loglike(output.prop)
      u<-runif(1)
      acc=min(1, exp(pi.new.log - pi.old.log), na.rm = T)
      acc.vec=c(acc.vec,acc)
      j=j+1
      
      if (u < acc) {
        output<-output.prop
      }
      output.mat<-rbind(output.mat, output)
      
      
      #################################################
      # Begin
      # update covariance matrix with adaptive MCMC
      #################################################
      if (j > 100) {
        if (j==101) {
          sigMat=cov(output.mat)
          thetaM=apply(output.mat, 2, mean)
        } else
        {
          tmp=update.cov(sigMat, i, thetaM, output)
          sigMat=tmp$sigMat
          thetaM=tmp$thetaM
        }
      }
      
      #################################################
      # End
      # update covariance matrix with adaptive MCMC
      #################################################
      
      
      
      ###############################################
      # Begin
      # update sigma using RM
      ###############################################
      if (m>n0) {
        sigma<-update.sigma(sigma2, acc, pstar, i)
        sigma2=sigma^2
        j=j+1
        sigma.vec<-c(sigma.vec, sigma)
        if ((j <= (iMax+n0)) && (Numbig<5 || Numsmall<5)) {
          Toobig<- (sigma > (3*sigma.start))
          Toosmall<-(sigma < (sigma.start/3))
          
          if (Toobig || Toosmall) {
            #restart the algorithm
            cat("restart the program at", i, "th iteration", "\n")
            sigma.restart<-c(sigma.restart, sigma)
            Numbig<- Numbig + Toobig
            Numsmall <- Numsmall + Toosmall
            i<-n0
            sigma.start<-sigma
          }
        } #end iMax
      }
      ###############################################
      # Begin
      # update sigma using RM
      ###############################################
      
      
    } #end niter
    #end of MCMC
    
    acc.vec <- acc.vec[ceiling(length(acc.vec)/2):length(acc.vec)]
    meanacc<-rep(NA, length(acc.vec))
    for (m in c(1:length(acc.vec))) {
      meanacc[m]=     mean(acc.vec[round(m/2) :m], na.rm = T)
    }
    # if(i == 1) {
    #   meanacc.record = meanacc
    # }
    meanacc.record[[i]] <- meanacc
    # plot(meanacc, type = "l", ylim = c(0, 0.7)); abline(h = 0.234)
    
    # for(m in 1:n.burnin) {
    #   output.prop<-as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
    #   pi.old.log <- loglike(output)
    #   pi.new.log <- loglike(output.prop)
    #   u<-runif(1)
    #   acc=min(1, exp(pi.new.log - pi.old.log))
    #   acc.vec=c(acc.vec,acc)
    #   if (u < acc) {
    #     output<-output.prop
    #   }
    # }
    output.prop<-as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
    pi.old.log <- loglike(output)
    # check
    if (is.na(pi.old.log)) {
      warning(paste("old log likelihood of patient ID =", newdata[[idVar]][1], "is NA!"))
    }
    pi.new.log <- loglike(output.prop)
    # check
    if (is.na(pi.new.log)) {
      warning(paste("new log likelihood of patient ID =", newdata[[idVar]][1], "is NA!"))
    }
    u<-runif(1)
    acc=min(1, exp(pi.new.log - pi.old.log))
    acc.vec=c(acc.vec,acc)
    #j=j+1
    
    if (u < acc) {
      output<-output.prop
    }
    output.mat<-rbind(output.mat, output)
  }
  
  # end ---------------------------------------------------------------------
  
  b <- output.mat[(nrow(output.mat)-iter+1):nrow(output.mat),]
  #
  # beta <- beta[(n.adapt+1):total,]
  # b <- b[(n.adapt+1):total,]
  # tau <- tau[(n.adapt+1):total,]
  # alpha <- alpha[(n.adapt+1):total,,,drop=F]
  # gamma <- gamma[(n.adapt+1):total,,,drop=F]
  # Sigma <- Sigma[(n.adapt+1):total,,]
  # gam_bh <- gam_bh[(n.adapt+1):total,,]
  
  # start calculating the risk ----------------------------------------------
  
  surv <- array(NA, iter)
  full.results <- list()
  overall.surv <- matrix(NA, iter, length(step))
  
  # denominator - survival part
  t.start <- step[1]
  XL.hs.now <- matrix(1,qp, ncol(pool.beta))
  ZL.hs.now <- matrix(1,qp, df.ns.longi + 1)
  XL.hs.now[,2:(df.ns.longi+1)] <- ZL.hs.now[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2,
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  }))
  XL.hs.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  
  XL.hs.back <- matrix(1, qp, ncol(pool.beta))
  ZL.hs.back <- matrix(1, qp, df.ns.longi + 1)
  XL.hs.back[,2:(df.ns.longi+1)] <- ZL.hs.back[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2 - 1,
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  }))
  XL.hs.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  
  T.hs <- splineDesign(knot.surv, (f.org + 1) * t/2, ord = 4L)
  
  hazard.survival <-array(NA, dim = c(iter, qp, 2))
  
  for (j in 1:iter) {
    hazard.survival[j,1:qp,1:2] <- exp(T.hs[1:qp,] %*% gam_bh[j,,1:2] + rep(gamma[j,,1:2] * X, each = qp) +
                                         (XL.hs.now[1:qp,] %*% beta[j,] + ZL.hs.now[1:qp,] %*% b[j,]) %*% alpha[j,1,1:2]  +
                                         (XL.hs.now[1:qp,] %*% beta[j,] + ZL.hs.now[1:qp,] %*% b[j,] -
                                            XL.hs.back[1:qp,] %*% beta[j,] - ZL.hs.back[1:qp,] %*% b[j,]) %*% alpha[j,2,1:2])
    surv[j] <-  exp(- t.start/2 * (rep(1,2) %*% (t(hazard.survival[j,1:qp,1:2]) %*% w)))
  }
  surv[surv < 1e-16] <- 1e-16
  
  for (i in 1:length(step)) {
    
    u <- step[i]
    t <- ifelse(i == 1, step[i], step[i-1])
    t.step <- f.org * (u-t)/2 + (u+t)/2
    T.h <- splineDesign(knot.surv, x = f.org * (u-t)/2 + (u+t)/2, ord = 4L)
    
    XL.h.now <- matrix(1, qp, ncol(pool.beta))
    ZL.h.now <- matrix(1, qp, df.ns.longi+1)
    XL.h.now[,2:(df.ns.longi+1)] <- ZL.h.now[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2,
                                                                       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    XL.h.back <- matrix(1, qp, ncol(pool.beta))
    ZL.h.back <- matrix(1, qp, df.ns.longi+1)
    XL.h.back[,2:(df.ns.longi+1)] <- ZL.h.back[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2 - 1,
                                                                         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    T.s <- array(1, dim = c(qp, ncol(pool.bs_gammas)/2, qp)) # second qp is the second integral
    XL.s.now <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.now <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.now[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:qp) {
      XL.s.now[j,2:(df.ns.longi+1),] <- ZL.s.now[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2,
                                                                               knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }
    
    XL.s.back <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.back <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.back[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:qp) {
      XL.s.back[j,2:(df.ns.longi+1),] <- ZL.s.back[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2 - 1,
                                                                                 knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }
    
    T.s.prep <- splineDesign(knot.surv, x = outer(f.org + 1, t.step/2), ord = 4L) # 1:qp - step[1]:f.org[1:qp]
    gk <- split(1:225, 1:qp)  # f.org group
    for (j in 1:qp) {
      T.s[,,j] <- T.s.prep[gk[[j]],]
    }
    
    # In addition (overall survival till different time points)
    XL.os.now <- matrix(1,qp, ncol(pool.beta))
    ZL.os.now <- matrix(1,qp, df.ns.longi + 1)
    XL.os.now[,2:(df.ns.longi+1)] <- ZL.os.now[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    XL.os.back <- matrix(1,qp, ncol(pool.beta))
    ZL.os.back <- matrix(1,qp, df.ns.longi + 1)
    XL.os.back[,2:(df.ns.longi+1)] <- ZL.os.back[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2 - 1,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    
    T.os <- splineDesign(knot.surv, (f.org + 1) * u/2, ord = 4L)
    
    
    survival <- matrix(NA, iter, qp)
    hazard <- matrix(NA, iter, qp)
    hazard.overall.survival <- array(NA, dim = c(iter, qp, 2))
    
    for (j in 1:iter) {
      for (m in 1:qp) {
        survival[j,m] <- exp(
          - t.step[m]/2 * sum(w %*% exp(t(T.s[m,,1:qp]) %*% gam_bh[j,,1:2]  + rep(gamma[j,,1:2] * X, each = qp) +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,]) %*% alpha[j,1,1:2] +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,] -
                                             t(XL.s.back[m,,1:qp]) %*% beta[j,] - t(ZL.s.back[m,,1:qp]) %*% b[j,]) %*% alpha[j,2,1:2]))
        )
        
      }
      hazard.overall.survival[j,1:qp,1:2] <- exp(T.os[1:qp,] %*% gam_bh[j,,1:2] + rep(gamma[j,,1:2] * X, each = qp) +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,]) %*% alpha[j,1,1:2] +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,] -
                                                      XL.os.back[1:qp,] %*% beta[j,] - ZL.os.back[1:qp,] %*% b[j,]) %*% alpha[j,2,1:2])
      hazard[j,1:qp] <- exp(T.h[1:qp,] %*% gam_bh[j,,1] +  gamma[j,,1] * X +
                              alpha[j,1,1] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,]) +
                              alpha[j,2,1] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,] -
                                                XL.h.back[1:qp,] %*% beta[j,] - ZL.h.back[1:qp,] %*% b[j,]))
    }
    risk.store[1:iter,i] <- (u-t)/2 * ((hazard[1:iter,] * survival[1:iter,]) %*% w) /
      c(surv[1:iter])
    
    overall.surv[1:iter,i] <- exp(- u/2 * (hazard.overall.survival[1:iter,,1] %*% w +
                                             hazard.overall.survival[1:iter,,2] %*% w))
    overall.surv[overall.surv[,i] < 1e-16,i] <- 1e-16
    # here we also record overall survival
    overall.surv.summary[i,2] <- mean(overall.surv[,i])
    overall.surv.summary[i,3:4] <- quantile(overall.surv[,i],
                                            probs= c(0.975, 0.025))
    
    risk[i,2] <- mean(rowSums(risk.store[,1:i, drop = F]))
    risk[i,3:4] <- quantile(rowSums(risk.store[,1:i, drop = F]),
                            probs= c(0.975, 0.025))
    full.results[[i]] <- list(rowSums(risk.store[,1:i, drop = F]))
  }
  
  #### Here are the treatment risk at the longest follow-up
  if (trt_part) {
    # matrices needed
    risk.store_trt <- rep(NA, iter)
    
    u <- max(Survtimes)
    t <- min(Survtimes)
    t.step <- f.org * (u-t)/2 + (u+t)/2
    T.h <- splineDesign(knot.surv, x = f.org * (u-t)/2 + (u+t)/2, ord = 4L)
    
    XL.h.now <- matrix(1, qp, ncol(pool.beta))
    ZL.h.now <- matrix(1, qp, df.ns.longi+1)
    XL.h.now[,2:(df.ns.longi+1)] <- ZL.h.now[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2,
                                                                       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    XL.h.back <- matrix(1, qp, ncol(pool.beta))
    ZL.h.back <- matrix(1, qp, df.ns.longi+1)
    XL.h.back[,2:(df.ns.longi+1)] <- ZL.h.back[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2 - 1,
                                                                         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    T.s <- array(1, dim = c(qp, ncol(pool.bs_gammas)/2, qp)) # second qp is the second integral
    XL.s.now <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.now <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.now[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:qp) {
      XL.s.now[j,2:(df.ns.longi+1),] <- ZL.s.now[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2,
                                                                               knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }
    
    XL.s.back <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.back <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.back[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]
    
    for (j in 1:qp) {
      XL.s.back[j,2:(df.ns.longi+1),] <- ZL.s.back[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2 - 1,
                                                                                 knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }
    
    T.s.prep <- splineDesign(knot.surv, x = outer(f.org + 1, t.step/2), ord = 4L) # 1:qp - step[1]:f.org[1:qp]
    gk <- split(1:225, 1:qp)  # f.org group
    for (j in 1:qp) {
      T.s[,,j] <- T.s.prep[gk[[j]],]
    }
    
    # In addition (overall survival till different time points)
    XL.os.now <- matrix(1,qp, ncol(pool.beta))
    ZL.os.now <- matrix(1,qp, df.ns.longi + 1)
    XL.os.now[,2:(df.ns.longi+1)] <- ZL.os.now[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    XL.os.back <- matrix(1,qp, ncol(pool.beta))
    ZL.os.back <- matrix(1,qp, df.ns.longi + 1)
    XL.os.back[,2:(df.ns.longi+1)] <- ZL.os.back[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2 - 1,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
    
    
    T.os <- splineDesign(knot.surv, (f.org + 1) * u/2, ord = 4L)
    
    
    survival <- matrix(NA, iter, qp)
    hazard <- matrix(NA, iter, qp)
    hazard.overall.survival <- array(NA, dim = c(iter, qp, 2))
    
    for (j in 1:iter) {
      for (m in 1:qp) {
        survival[j,m] <- exp(
          - t.step[m]/2 * sum(w %*% exp(t(T.s[m,,1:qp]) %*% gam_bh[j,,1:2]  + rep(gamma[j,,1:2] * X, each = qp) +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,]) %*% alpha[j,1,1:2] +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,] -
                                             t(XL.s.back[m,,1:qp]) %*% beta[j,] - t(ZL.s.back[m,,1:qp]) %*% b[j,]) %*% alpha[j,2,1:2]))
        )
        
      }
      hazard.overall.survival[j,1:qp,1:2] <- exp(T.os[1:qp,] %*% gam_bh[j,,1:2] + rep(gamma[j,,1:2] * X, each = qp) +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,]) %*% alpha[j,1,1:2] +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,] -
                                                      XL.os.back[1:qp,] %*% beta[j,] - ZL.os.back[1:qp,] %*% b[j,]) %*% alpha[j,2,1:2])
      hazard[j,1:qp] <- exp(T.h[1:qp,] %*% gam_bh[j,,2] +  gamma[j,,2] * X +
                              alpha[j,1,2] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,]) +
                              alpha[j,2,2] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,] -
                                                XL.h.back[1:qp,] %*% beta[j,] - ZL.h.back[1:qp,] %*% b[j,]))
    }
    risk.store_trt <- (u-t)/2 * ((hazard[1:iter,] * survival[1:iter,]) %*% w) /
      c(surv[1:iter])
    
    overall.surv <- exp(- u/2 * (hazard.overall.survival[1:iter,,1] %*% w +
                                   hazard.overall.survival[1:iter,,2] %*% w))
    overall.surv[overall.surv < 1e-16] <- 1e-16
    risk_trt <- data.frame(t.hor = max(Survtimes),
                       mean = NA,
                       upper = NA,
                       lower = NA)
    risk_trt[,2] <- mean(risk.store_trt)
    risk_trt[,3:4] <- quantile(risk.store_trt,
                            probs= c(0.975, 0.025))
  }
  
  longi[,2] <- sapply(seq(0, max(test$TimeSince_Dx), length.out = 50), function(n) {
    mean(sapply(1:iter, function(l) {
      c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), test.id$DxAge) %*% beta[l,] +
        c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)])) %*% b[l,]
    }))
  })
  
  if (trt_part) {
    inf_mean <- mean(rowSums(risk.store)/(rowSums(risk.store) + risk.store_trt))
    inf_ci <- quantile(rowSums(risk.store)/(rowSums(risk.store) + risk.store_trt), 
                       probs= c(0.975, 0.025))
    return(list(risk_visits = rbind(risk[risk$t.hor %in% c(t_start, t_visits),],
                                    c(Inf,inf_mean, inf_ci[1], inf_ci[2])),
                risk = risk,
                risk_trt = risk_trt,
                longi = longi,
                acc = meanacc.record,
                surv = surv,
                parameters = list(beta = beta,
                                  b = b,
                                  gam_bh = gam_bh,
                                  gamma = gamma,
                                  alpha = alpha,
                                  Sigma = Sigma,
                                  tau = tau),
                inits = list(sigma = sigma, sigMat = sigMat),
                overall.surv = overall.surv,
                overall.surv.summary = overall.surv.summary,
                full.results = full.results))
  } else {
    return(list(risk_visits = risk[risk$t.hor %in% c(t_start, t_visits),],
                risk = risk,
                longi = longi,
                acc = meanacc.record,
                surv = surv,
                parameters = list(beta = beta,
                                  b = b,
                                  gam_bh = gam_bh,
                                  gamma = gamma,
                                  alpha = alpha,
                                  Sigma = Sigma,
                                  tau = tau),
                inits = list(sigma = sigma, sigMat = sigMat),
                overall.surv = overall.surv,
                overall.surv.summary = overall.surv.summary,
                full.results = full.results))
  }
  
}
