##----------------------------individual bias correction----------------------------##
# try this again because least square is minimum variance estimator...?
# individual bias correction by regressing Y on the mean of each point, 
# each Y is a realization of dis
# this is by definition of bias

# new simple bias correction (mean error)
## normal dis => means = median
mean_err<-function(Y, means){
  emp_err <- means-Y
  err_means <- colMeans(emp_err)
  return(err_means)
}
##------------------------------Pooling Functions -----------------------------------##
# individual recalibration 1
wrapper_indBLP <- function(Y, means, sds, pars, eval = TRUE,varPIT = FALSE){
  ab     <- exp(pars[1:2])
  g      <- dnorm(Y, mean = means, sd = sds) 
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds), 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds), shape1 = ab[1], shape2 = ab[2])
  if(varPIT) return(var(PIT))
  if(eval){
    return(sum(ll))
  }else{
    return(ll)
  }
}

# ensemble calibration
wrapper_TLP <- function(Y, means, sds, weights, eval = TRUE, 
                        plot = FALSE){
  weights<-exp(weights)
  weights <- weights/sum(weights)
  ll <- log(dnorm(Y, mean = means, sd = sds) %*% weights)
  PIT <- pnorm(Y, mean = means, sd = sds) %*% weights
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "TLP")
    abline( h = 1, lty = 2)
  }
  if(eval){
    return(sum(ll))
  } 
  else{
    return(ll)
  }
}


wrapper_BLP <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE){
  
  ab     <- exp(pars[1:2])
  weights <- exp(pars[-c(1,2)])
  weights <- weights/sum(weights)
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "BLP")
    abline( h = 1, lty = 2)
  }  
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}


wrapper_nBLP <- function(Y, means, sds, pars, eval = TRUE, 
                         plot = FALSE){
  
  abc     <- exp(pars[1:3])
  weights <- exp(pars[-c(1:3)])
  
  weights <- weights/sum(weights)
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "Non-central BLP")
    abline( h = 1, lty = 2)
  }  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}


wrapper_cBLP <- function(Y, means, sds, pars, cab, eval = TRUE, 
                         plot = FALSE){
  comp_bpars<-cab
  # normal probs
  comp_nprobs<-cbind(pnorm(Y, mean = means[,1], sd = sds[,1]),
                     pnorm(Y, mean = means[,2], sd = sds[,2]),
                     pnorm(Y, mean = means[,3], sd = sds[,3]))
  # get beta probs
  comp_bprobs<-cbind(pbeta(pnorm(Y, mean = means[,1], sd = sds[,1]), 
                           shape1 = comp_bpars[1,1], shape2 = comp_bpars[1,2]),
                     pbeta(pnorm(Y, mean = means[,2], sd = sds[,2]), 
                           shape1 = comp_bpars[2,1], shape2 = comp_bpars[2,2]),
                     pbeta(pnorm(Y, mean = means[,3], sd = sds[,3]), 
                           shape1 = comp_bpars[3,1], shape2 = comp_bpars[3,2]))
  ab     <- exp(pars[1:2])
  weights <- exp(pars[-c(1,2)])
  
  weights <- weights/sum(weights)
  f1 <- as.matrix(dnorm(Y, mean = means, sd = sds))
  f2<- as.matrix(dbeta(comp_nprobs,shape1 = comp_bpars[,1], shape2 = comp_bpars[,2]))
  g      <- (f1*f2) %*% weights
  g_beta <- dbeta(comp_bprobs %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(comp_bprobs %*% weights, shape1 = ab[1], shape2 = ab[2])
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "cBLP")
    abline( h = 1, lty = 2)
  }  
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

# wrapper_cBLP2 <- function(Y, means, sds, pars, cab, var, eval = TRUE, 
#                          plot = FALSE){
#   comp_bpars<-cab
#   # normal probs
#   comp_nprobs<-cbind(pnorm(Y, mean = means[,1], sd = sds[,1]),
#                      pnorm(Y, mean = means[,2], sd = sds[,2]),
#                      pnorm(Y, mean = means[,3], sd = sds[,3]))
#     
#   # get beta probs
#   for(i in 1:3){
#     if(var[i]<(1/12)){
#       p1 <- 
#     }
#     comp_nprobs<-cbind(pnorm(Y, mean = means[,1], sd = sds[,1]),
#                        pnorm(Y, mean = means[,2], sd = sds[,2]),
#                        pnorm(Y, mean = means[,3], sd = sds[,3]))
#     
#   }
#   comp_nprobs <- cbind(p1,p2,p3)
#   
#   comp_bprobs<-cbind(pbeta(pnorm(Y, mean = means[,1], sd = sds[,1]), 
#                            shape1 = comp_bpars[1,1], shape2 = comp_bpars[1,2]),
#                      pbeta(pnorm(Y, mean = means[,2], sd = sds[,2]), 
#                            shape1 = comp_bpars[2,1], shape2 = comp_bpars[2,2]),
#                      pbeta(pnorm(Y, mean = means[,3], sd = sds[,3]), 
#                            shape1 = comp_bpars[3,1], shape2 = comp_bpars[3,2]))
#   
#   ab     <- exp(pars[1:2])
#   weights <- exp(pars[-c(1,2)])
#   
#   weights <- weights/sum(weights)
#   f1 <- as.matrix(dnorm(Y, mean = means, sd = sds))
#   f2<- as.matrix(dbeta(comp_nprobs,shape1 = comp_bpars[,1], shape2 = comp_bpars[,2]))
#   g      <- (f1*f2) %*% weights
#   g_beta <- dbeta(comp_bprobs %*% weights, 
#                   shape1 = ab[1], shape2 = ab[2])
#   
#   ll     <- log(g) + log(g_beta)
#   
#   PIT <- pbeta(comp_bprobs %*% weights, shape1 = ab[1], shape2 = ab[2])
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.8), main = "cBLP2")
#     abline( h = 1, lty = 2)
#   }  
#   
#   if(eval){
#     return(sum(ll))
#   } else{
#     return(ll)
#   }
# }

# while function 
# individual ones
get_indBLPpars<-function(means,sds,modelnum){
  ab          <- log(c(1,1))
  #  while(TRUE){
  est <- optim(
    par     = ab,
    fn      = wrapper_indBLP,
    Y       = Y[ind_est],
    means   = means[ind_est,modelnum],
    sds   = sds[ind_est,modelnum],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE, maxit = 200)
  )
  ab          <- exp(est$par)
  #   if(est$par>0) break
  #  }
  return(ab)
}

# ensemble
get_TLPpars<-function(means,sds){
  weights_TLP <- log(c(0.3, 0.3, 0.4))
  #  while(TRUE){
  est <- optim(
    par     = weights_TLP,
    fn      = wrapper_TLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE, maxit = 200)
  )
  weights_TLP <- exp(est$par)/sum(exp(est$par))
  #    if(abs(sum(est$par) - 1) < 0.01) break
  #  }
  return(weights_TLP)
}


get_BLPpars<-function(means,sds){
  weights_BLP <- log(c(0.3, 0.3, 0.4))
  ab          <- log(c(1,1))
  pars_BLP    <- c(ab, weights_BLP)
  
  # Estimate weights for BLP
  #  while(TRUE){
  
  est <- optim(
    par     = pars_BLP,
    fn      = wrapper_BLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE, maxit = 200)
  )
  
  ab          <- exp(est$par[1:2])
  weights_BLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
  pars_BLP    <- c(ab, weights_BLP)
  #    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
  #  }
  return(pars_BLP)
}

get_nBLPpars<-function(means,sds){
  weights_nBLP <- log(c(0.3, 0.3, 0.4))
  abc          <- log(c(1,1,0.01))
  pars_nBLP    <- c(abc, weights_nBLP)
  #  while(TRUE){
  est <- optim(
    par     = pars_nBLP,
    fn      = wrapper_nBLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    method  = "BFGS",
    control = list(fnscale = -1,trace = TRUE, maxit = 200)
    
  )
  
  abc          <- exp(est$par[1:3])
  weights_nBLP <- exp(est$par[-c(1:3)])/sum(exp(est$par[-c(1:3)]))
  pars_nBLP    <- c(abc, weights_nBLP)
  
  #    if(abs(sum(est$par[-c(1:3)]) - 1) < 0.01) break
  #  }
  return(pars_nBLP)
}

get_cBLPpars<-function(means,sds,cab){
  weights_CBLP <- log(c(0.3, 0.3, 0.4))
  ab          <- log(c(1,1))
  pars_CBLP    <- c(ab, weights_CBLP)
  #  while(TRUE){
  
  est <- optim(
    par     = pars_CBLP,
    fn      = wrapper_cBLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    cab   = cab,
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE, maxit = 200)
  )
  ab          <- exp(est$par[1:2])
  weights_CBLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
  pars_CBLP    <- c(ab, weights_CBLP)
  #    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
  #  }
  return(pars_CBLP)
}

##---------------------- PIT and log likelihood functions--------------------------------##
PIT_normal <- function(Y, mu, sd, var = FALSE, plot = FALSE,name){
  PIT <- pnorm(Y, mean = mu, sd = sd)
  
  if(var) return(var(PIT))
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 2), main = paste0("Density Forecast",name))
    abline( h = 1, lty = 2)
  }}


ll_normal         <- function(Y, mu, sd){
  ll <- log(dnorm(Y, mean = mu, sd = sd))
  return(ll)
}

wrapper_indBLP_nonc <- function(Y, means, sds, pars, eval = TRUE,plot = FALSE,
                                varPIT = FALSE, name){
  
  ab     <- pars[1:2]
  g      <- dnorm(Y, mean = means, sd = sds) 
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds), 
                  shape1 = ab[1], shape2 = ab[2])
  ll     <- log(g) + log(g_beta)
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds), shape1 = ab[1], shape2 = ab[2])
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = paste0("BT Density Forecast",name))
    abline( h = 1, lty = 2)
  }
  if(varPIT) return(var(PIT))
  if(eval){
    return(sum(ll))
  }else{
    return(ll)
  }
}

wrapper_TLP_nonc <- function(Y, means, sds, weights, eval = TRUE, plot = FALSE, varPIT=FALSE){
  ll <- log(dnorm(Y, mean = means, sd = sds) %*% weights)
  PIT <- pnorm(Y, mean = means, sd = sds) %*% weights
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "TLP")
    abline( h = 1, lty = 2)
  }
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } 
  else{
    return(ll)
  }
}


wrapper_BLP_nonc  <- function(Y, means, sds, pars, eval = TRUE, plot = FALSE,varPIT=FALSE){
  
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]
  
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "BLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

# these two wrappers below are for plots only
wrapper_nBLP_nonc  <- function(Y, means, sds, pars, eval = TRUE, plot = FALSE,varPIT=FALSE){
  
  abc     <- pars[1:3]
  weights <- pars[-c(1:3)]
  g      <- dnorm(Y, mean = means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds) %*% weights, 
                  shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = means, sd = sds) %*% weights, shape1 = abc[1], shape2 = abc[2],ncp=abc[3])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "Non-central BLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

wrapper_cBLP_nonc  <- function(Y, means, sds, pars, cab, eval = TRUE, 
                               plot = FALSE, varPIT = FALSE){
  comp_bpars<-cab
  # normal probs
  comp_nprobs<-cbind(pnorm(Y, mean = means[,1], sd = sds[,1]),
                     pnorm(Y, mean = means[,2], sd = sds[,2]),
                     pnorm(Y, mean = means[,3], sd = sds[,3]))
  # get beta probs
  comp_bprobs<-cbind(pbeta(pnorm(Y, mean = means[,1], sd = sds[,1]), 
                           shape1 = comp_bpars[1,1], shape2 = comp_bpars[1,2]),
                     pbeta(pnorm(Y, mean = means[,2], sd = sds[,2]), 
                           shape1 = comp_bpars[2,1], shape2 = comp_bpars[2,2]),
                     pbeta(pnorm(Y, mean = means[,3], sd = sds[,3]), 
                           shape1 = comp_bpars[3,1], shape2 = comp_bpars[3,2]))
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]
  f1 <- as.matrix(dnorm(Y, mean = means, sd = sds))
  f2<- as.matrix(dbeta(comp_nprobs,shape1 = comp_bpars[,1], shape2 = comp_bpars[,2]))
  g      <- (f1*f2) %*% weights
  g_beta <- dbeta(comp_bprobs %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(comp_bprobs %*% weights, shape1 = ab[1], shape2 = ab[2])
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "cBLP")
    abline( h = 1, lty = 2)
  }  
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

##---------------------- wrappers --------------------------------##

make_table <- function(forecast_names,TLP_pars,BLP_pars,nBLP_pars,cBLP_pars){
  table <- matrix(NA, nrow = 4, ncol = 6, 
                  dimnames = list(c(forecast_names), c("w1", "w2", "w3", "alpha", "beta","ncp")))
  table[1,] <- c(TLP_pars, rep(NA, 3))
  table[2,] <- c(BLP_pars[-c(1:2)], BLP_pars[1:2],NA)
  table[3,] <- c(nBLP_pars[-c(1:3)], nBLP_pars[1:3])
  table[4,] <- c(cBLP_pars[-c(1:2)], cBLP_pars[1:2],NA)
  table<-round(table, 3)
  return(table)
}

make_ls <- function(forecast_names,Y,means,sds,TLP_pars,BLP_pars,nBLP_pars,cBLP_pars1,cBLP_pars2){
  mls <- matrix(NA, nrow = 7, ncol = 2,
                 dimnames = list(c(c("f1","f2","f3"),forecast_names),
                                 c("Training", "Test")))
  for(j in 1:2){
    ind <- if(j == 1){ind_est}else{!ind_est}
    for(i in 1:3) mls[i, j]  <- mean(ll_normal(Y = Y[ind], mu = means[ind, i],
                                                sd = sds[ind, i]))
    mls[4, j] <- mean(wrapper_TLP_nonc(Y = Y[ind], means = means[ind,],
                                        sds = sds[ind,], weights = TLP_pars, eval = FALSE))
    mls[5, j] <- mean(wrapper_BLP_nonc(Y = Y[ind], means = means[ind,],
                                        sds = sds[ind,], pars = BLP_pars, eval = FALSE))
    mls[6, j] <- mean(wrapper_nBLP_nonc(Y = Y[ind], means = means[ind,],
                                         sds = sds[ind,], pars = nBLP_pars, eval = FALSE))
    mls[7, j] <- mean(wrapper_cBLP_nonc(Y = Y[ind], means = means[ind,], sds = sds[ind,],
                                         pars = cBLP_pars1,cab=cBLP_pars2, eval = FALSE))
  }
  mls<-round(mls, 3)
  return(mls)
}

make_var <- function(forecast_names,Y,means,sds,TLP_pars,BLP_pars,nBLP_pars,cBLP_pars1,cBLP_pars2){
  var <- matrix(NA, nrow = 7, ncol = 2,
                 dimnames = list(c(c("f1","f2","f3"),forecast_names),
                                 c("Training", "Test")))
  for(j in 1:2){
    ind <- if(j == 1){ind_est}else{!ind_est}
    for(i in 1:3) var[i, j]  <- mean(PIT_normal(Y = Y[ind], mu = means[ind, i],
                                                 sd = sds[ind, i],var=TRUE))
    var[4, j] <- mean(wrapper_TLP_nonc(Y = Y[ind], means = means[ind,],
                                        sds = sds[ind,], weights = TLP_pars, eval = FALSE,
                                        varPIT=TRUE))
    var[5, j] <- mean(wrapper_BLP_nonc(Y = Y[ind], means = means[ind,],
                                        sds = sds[ind,], pars = BLP_pars, eval = FALSE,
                                        varPIT=TRUE))
    var[6, j] <- mean(wrapper_nBLP_nonc(Y = Y[ind], means = means[ind,],
                                         sds = sds[ind,], pars = nBLP_pars, eval = FALSE,
                                         varPIT=TRUE))
    var[7, j] <- mean(wrapper_cBLP_nonc(Y = Y[ind], means = means[ind,], sds = sds[ind,],
                                         pars = cBLP_pars1,cab=cBLP_pars2, eval = FALSE,
                                         varPIT=TRUE))
  }
  var<-round(var, 3)
  return(var)
}

plot_PITs <- function(Y,means,sds,TLP_pars,BLP_pars,nBLP_pars,cBLP_pars1,cBLP_pars2){
  for(i in 1:3) PIT_normal(Y = Y[!ind_est], mu = means[!ind_est ,i], 
                           sd = sds[!ind_est ,i], plot = TRUE, name = i)
  for(i in 1:3) wrapper_indBLP_nonc(Y = Y[!ind_est], means = means[!ind_est ,i], 
                                    sds = sds[!ind_est ,i], pars = cBLP_pars2[i,],
                                    plot = TRUE, name = i) 
  wrapper_TLP_nonc(Y = Y[!ind_est], means = means[!ind_est, ],
                   sds = sds[!ind_est, ], weights = TLP_pars, plot = TRUE)
  wrapper_BLP_nonc(Y = Y[!ind_est], means = means[!ind_est, ],
                   sds = sds[!ind_est, ], pars = BLP_pars, plot = TRUE)
  wrapper_nBLP_nonc(Y = Y[!ind_est], means = means[!ind_est, ],
                    sds = sds[!ind_est, ], pars = nBLP_pars, plot = TRUE)
  wrapper_cBLP_nonc(Y = Y[!ind_est], means = means[!ind_est, ],
                    sds = sds[!ind_est, ], pars = cBLP_pars1, cab=cBLP_pars2, plot = TRUE)
}
# wrapper_SLP_nonc <- function(Y, means, sds, pars, eval = TRUE, 
#                              plot = FALSE,varPIT=FALSE){
#   
#   con     <- pars[1]
#   weights <- pars[-1]
#   # normal - mean and med the same
#   means_SLP <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)/con
#   
#   ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
#   
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "SLP")
#     abline( h = 1, lty = 2)
#   }
#   if(varPIT) return(var(PIT))
#   
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }

# wrapper_cSLP_nonc <- function(Y, means, sds, pars, ind_con, eval = TRUE, 
#                          plot = FALSE,varPIT=FALSE){
#   
#   con     <- pars[1]
#   weights <- pars[-1]
#   # normal - mean and med the same
#   means_SLP <- matrix(rep(Y, ncol(means)), ncol = ncol(means))/con
#   ll <- log(1/con * dnorm(means_SLP,mean = 0,sd=sds/as.numeric(ind_con)) %*% weights)
#   
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds/as.numeric(ind_con)) %*% weights
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "cSLP")
#     abline( h = 1, lty = 2)
#   }
#   if(varPIT) return(var(PIT))
#   
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }

#### add back SLP

# wrapper_SLP <- function(Y, means, sds, pars, eval = TRUE, 
#                         plot = FALSE){
#   
#   con     <- exp(pars[1])
#   weights <- exp(pars[-1])
#   weights <- weights/sum(weights)
#   # normal - mean and med the same
#   means_SLP <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)/con
#   
#   ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
#   
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "SLP")
#     abline( h = 1, lty = 2)
#   }
#   
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }
# 
# get_SLPpars<-function(means,sds){
#   weights_SLP <- log(c(0.3, 0.3, 0.4))
#   c         <- log(1.000000000000001)
#   pars_SLP    <- c(c, weights_SLP)
#   est <- optim(
#     par     = pars_SLP,
#     fn      = wrapper_SLP,
#     Y       = Y[ind_est],
#     means   = means[ind_est, ],
#     sds   = sds[ind_est, ],
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   c          <- exp(est$par[1])
#   weights_SLP <- exp(est$par[-1])/sum(exp(est$par[-1]))
#   pars_SLP    <- c(c, weights_SLP)
#   return(pars_SLP)
# }
# 
# ## check
# 
# wrapper_cSLP <- function(Y, means, sds, pars, ind_con, eval = TRUE, 
#                         plot = FALSE){
#   
#   con     <- exp(pars[1])
#   weights <- exp(pars[-1])
#   weights <- weights/sum(weights)
#   # normal - mean and med the same
#   means_SLP <- matrix(rep(Y, ncol(means)), ncol = ncol(means)) /con
#   ll <- log(1/con * dnorm(means_SLP,mean = 0,sd=sds/as.numeric(ind_con)) %*% weights)
#   
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds/as.numeric(ind_con)) %*% weights
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "cSLP")
#     abline( h = 1, lty = 2)
#   }
#   
# 
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }
# 
# get_cSLPpars<-function(means,sds,ind_con){
#   weights_SLP <- log(c(0.3, 0.3, 0.4))
#   c         <- log(1.000000000000001)
#   pars_SLP    <- c(c, weights_SLP)
#   est <- optim(
#     par     = pars_SLP,
#     fn      = wrapper_cSLP,
#     Y       = Y[ind_est],
#     means   = means[ind_est, ],
#     sds   = sds[ind_est, ],
#     ind_con = ind_con,
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   c          <- exp(est$par[1])
#   weights_SLP <- exp(est$par[-1])/sum(exp(est$par[-1]))
#   pars_SLP    <- c(c, weights_SLP)
#   return(pars_SLP)
# }
# 
# #### individual SLP
# 
# wrapper_indSLP <- function(Y, means, sds, pars, eval = TRUE, 
#                         plot = FALSE){
#   
#   con     <- exp(pars)
#   means_SLP <- (Y - means)/con
#   ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds))
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds)
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "SLP")
#     abline( h = 1, lty = 2)
#   }
#   
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }
# 
# get_indSLPpars<-function(means,sds,modelnum){
#   pars_indSLP <- log(1.000000000000001)
#   est <- optim(
#     par     = pars_indSLP,
#     fn      = wrapper_indSLP,
#     Y       = Y[ind_est],
#     means   = means[ind_est,modelnum],
#     sds   = sds[ind_est,modelnum],
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   pars_indSLP <- exp(est$par)
#   return(pars_indSLP)
# }
# these two wrappers below are for plots only
# wrapper_bcBLP <- function(Y, means, sds, pars, errs, eval = TRUE, 
#                           plot = FALSE){
#   ab     <- exp(pars[1:2])
#   weights <- exp(pars[-c(1,2)])
#   weights <- weights/sum(weights)
#   corrected_means <- means + errs
#   g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
#   g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                   shape1 = ab[1], shape2 = ab[2])
#   
#   ll     <- log(g) + log(g_beta)
#   
#   PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.8), main = "bcBLP")
#     abline( h = 1, lty = 2)
#   }
#   if(eval){
#     return(sum(ll))
#   } 
#   else{
#     return(ll)
#   }
# }
# wrapper_bcnBLP <- function(Y, means, sds, pars, errs, eval = TRUE, 
#                            plot = FALSE){
#   ab     <- exp(pars[1:3])
#   weights <- exp(pars[-c(1:3)])
#   weights <- weights/sum(weights)
#   corrected_means <- means + errs
#   g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
#   g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                   shape1 = ab[1], shape2 = ab[2], ncp = ab[3])
#   
#   ll     <- log(g) + log(g_beta)
#   
#   PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                shape1 = ab[1], shape2 = ab[2],ncp=ab[3])
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.8), main = "bcnBLP")
#     abline( h = 1, lty = 2)
#   }
#   if(eval){
#     return(sum(ll))
#   } 
#   else{
#     return(ll)
#   }
# }
# bcBLP
# get_bcBLPpars<-function(means,sds, errs){
#   weights_bcBLP <- log(c(0.3, 0.3, 0.4))
#   ab          <- log(c(1,1))
#   pars_bcBLP    <- c(ab, weights_bcBLP)
#   
#   # Estimate weights for BLP
#   #  while(TRUE){
#   
#   est <- optim(
#     par     = pars_bcBLP,
#     fn      = wrapper_bcBLP,
#     Y       = Y[ind_est],
#     means   = means[ind_est, ],
#     sds   = sds[ind_est, ],
#     errs = errs,
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   
#   ab          <- exp(est$par[1:2])
#   weights_bcBLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
#   pars_bcBLP    <- c(ab, weights_bcBLP)
#   #    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
#   #  }
#   return(pars_bcBLP)
# }
# get_bcnBLPpars<-function(means,sds, errs){
#   weights_bcnBLP <- log(c(0.3, 0.3, 0.4))
#   ab          <- log(c(1,1,0.01))
#   pars_bcnBLP    <- c(ab, weights_bcnBLP)
#   
#   # Estimate weights for BLP
#   #  while(TRUE){
#   
#   est <- optim(
#     par     = pars_bcnBLP,
#     fn      = wrapper_bcnBLP,
#     Y       = Y[ind_est],
#     means   = means[ind_est, ],
#     sds   = sds[ind_est, ],
#     errs = errs,
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   
#   ab          <- exp(est$par[1:3])
#   weights_bcnBLP <- exp(est$par[-c(1:3)])/sum(exp(est$par[-c(1:3)]))
#   pars_bcnBLP    <- c(ab, weights_bcnBLP)
#   return(pars_bcnBLP)
# }
# wrapper_bcnBLP_nonc  <- function(Y, means, sds, pars, errs, eval = TRUE, plot = FALSE,varPIT=FALSE){
#   
#   ab     <- pars[1:3]
#   weights <- pars[-c(1:3)]
#   corrected_means <- means + errs
#   g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
#   g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                   shape1 = ab[1], shape2 = ab[2], ncp = ab[3])
#   
#   ll     <- log(g) + log(g_beta)
#   
#   PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                shape1 = ab[1], shape2 = ab[2], ncp = ab[3])
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.8), main = "bcnBLP")
#     abline( h = 1, lty = 2)
#   } 
#   if(varPIT) return(var(PIT))
#   
#   if(eval){
#     return(sum(ll))
#   } else{
#     return(ll)
#   }
# }
# wrapper_bcBLP_nonc  <- function(Y, means, sds, pars, errs, eval = TRUE, plot = FALSE,varPIT=FALSE){
#   
#   ab     <- pars[1:2]
#   weights <- pars[-c(1,2)]
#   corrected_means <- means + errs
#   g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
#   g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
#                   shape1 = ab[1], shape2 = ab[2])
#   
#   ll     <- log(g) + log(g_beta)
#   
#   PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.8), main = "bcBLP")
#     abline( h = 1, lty = 2)
#   } 
#   if(varPIT) return(var(PIT))
#   
#   if(eval){
#     return(sum(ll))
#   } else{
#     return(ll)
#   }
# }
# wrapper_SLP_nonc <- function(Y, means, sds, pars, eval = TRUE, 
#                              plot = FALSE,varPIT=FALSE){
#   
#   con     <- pars[1]
#   weights <- pars[-1]
#   # normal - mean and med the same
#   means_SLP <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)/con
#   
#   ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
#   
#   PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
#   
#   if(plot){
#     hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
#          ylim = c(0, 1.3), main = "SLP")
#     abline( h = 1, lty = 2)
#   }
#   if(varPIT) return(var(PIT))
#   
#   if(eval){
#     return(sum(ll))
#   }else{
#     return(ll)
#   }
# }
# get_ind_cbcBLPpars<-function(means,sds,modelnum){
#   ab          <- log(c(1,1))
#   #  while(TRUE){
#   est <- optim(
#     par     = ab,
#     fn      = wrapper_indBLP,
#     Y       = Y[ind_est],
#     means   = means[,modelnum],
#     sds   = sds[ind_est,modelnum],
#     method  = "BFGS",
#     control = list(fnscale = -1, trace = TRUE, maxit = 200)
#   )
#   ab          <- exp(est$par)
#   #   if(est$par>0) break
#   #  }
#   return(ab)
# }