##----------------------------individual bias correction----------------------------##
# try this again because least square is minimum variance estimator...?
# individual bias correction by regressing Y on the mean of each point, 
# each Y is a realization of dis
# this is by definition of bias
# 
# # to get fitted (biased corrected values)
bc_slr<-function(Y, means){
  Yind<-Y[ind_est]
  means_ind<-means[ind_est, ]
  g1<-lm(Yind~means_ind[,1])$fitted.values
  g2<-lm(Yind~means_ind[,2])$fitted.values
  g3<-lm(Yind~means_ind[,3])$fitted.values
  fittedVal_g<-as.matrix(cbind(g1,g2,g3),ncol=3)
  return(fittedVal_g)
}
mean_err<-function(Y, means){
  emp_err <- means-Y
  err_means <- colMeans(emp_err)
  return(err_means)
}
bc_ests<-function(Y, means){
  Yind<-Y[ind_est]
  means_ind<-means[ind_est, ]
  g1e<-lm(Yind~means_ind[,1])$coefficients
  g2e<-lm(Yind~means_ind[,2])$coefficients
  g3e<-lm(Yind~means_ind[,3])$coefficients
  coef_g<-as.matrix(rbind(g1e,g2e,g3e),ncol=2)
  return(coef_g)
}

# new simple bias correction (mean error)
## normal dis => means = median
mean_err<-function(Y, means){
  emp_err <- means-Y
  err_means <- colMeans(emp_err)
  return(err_means)
}
##------------------------------Pooling Functions -----------------------------------##

# individual recalibration
wrapper_indBLP <- function(Y, means, sds, pars, eval = TRUE){
  
  ab     <- exp(pars[1:2])
  
  g      <- dnorm(Y, mean = means, sd = sds) 
  g_beta <- dbeta(pnorm(Y, mean = means, sd = sds), 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
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
# these two wrappers below are for plots only
wrapper_bcBLP <- function(Y, means, sds, pars, errs, eval = TRUE, 
                        plot = FALSE){
  ab     <- exp(pars[1:2])
  weights <- exp(pars[-c(1,2)])
  weights <- weights/sum(weights)
  corrected_means <- means + errs
  g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "bcBLP")
    abline( h = 1, lty = 2)
  }
  if(eval){
    return(sum(ll))
  } 
  else{
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

wrapper_cnBLP <- function(Y, means, sds, pars, eval = TRUE, plot = FALSE){
  
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
         ylim = c(0, 1.8), main = "cNBLP")
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
      control = list(fnscale = -1, trace = TRUE)
    )
  ab          <- exp(est$par)
#   if(est$par>0) break
#  }
  return(ab)
}

get_indNBLPpars<-function(means,sds,modelnum){
  abc          <- log(c(1,1))
  est <- optim(
    par     = abc,
    fn      = wrapper_indNBLP,
    Y       = Y[ind_est],
    means   = means[ind_est,modelnum],
    sds   = sds[ind_est,modelnum],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE)
  )
  abc          <- exp(est$par)
  return(abc)
}

get_ind_cbcBLPpars<-function(means,sds,modelnum){
  ab          <- log(c(1,1))
  #  while(TRUE){
  est <- optim(
    par     = ab,
    fn      = wrapper_indBLP,
    Y       = Y[ind_est],
    means   = means[,modelnum],
    sds   = sds[ind_est,modelnum],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE)
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
      control = list(fnscale = -1, trace = TRUE)
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
      control = list(fnscale = -1, trace = TRUE)
    )
    
    ab          <- exp(est$par[1:2])
    weights_BLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
    pars_BLP    <- c(ab, weights_BLP)
#    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
#  }
  return(pars_BLP)
}
# bcTLP
get_bcTLPpars<-function(means,sds){
  weights_TLP <- log(c(0.3, 0.3, 0.4))
  #  while(TRUE){
  est <- optim(
    par     = weights_TLP,
    fn      = wrapper_TLP,
    Y       = Y[ind_est, ],
    means   = means,
    sds   = sds[ind_est, ],
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE)
  )
 weights_TLP <- exp(est$par)/sum(exp(est$par))

  #    if(abs(sum(est$par) - 1) < 0.01) break
  #  }
  return(weights_TLP)
}

# bcBLP
get_bcBLPpars<-function(means,sds, errs){
  weights_bcBLP <- log(c(0.3, 0.3, 0.4))
  ab          <- log(c(1,1))
  pars_bcBLP    <- c(ab, weights_bcBLP)
  
  # Estimate weights for BLP
  #  while(TRUE){
  
  est <- optim(
    par     = pars_bcBLP,
    fn      = wrapper_bcBLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    errs = errs,
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE)
  )
  
  ab          <- exp(est$par[1:2])
  weights_bcBLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
  pars_bcBLP    <- c(ab, weights_bcBLP)
  #    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
  #  }
  return(pars_bcBLP)
}

get_nBLPpars<-function(means,sds){
  weights_nBLP <- log(c(0.3, 0.3, 0.4))
  abc          <- log(c(1,1,0.1))
  pars_nBLP    <- c(abc, weights_nBLP)
#  while(TRUE){
    est <- optim(
      par     = pars_nBLP,
      fn      = wrapper_nBLP,
      Y       = Y[ind_est],
      means   = means[ind_est, ],
      sds   = sds[ind_est, ],
      method  = "BFGS",
      control = list(fnscale = -1,trace = TRUE)

    )

    abc          <- exp(est$par[1:3])
    weights_nBLP <- exp(est$par[-c(1:3)])/sum(exp(est$par[-c(1:3)]))
    pars_nBLP    <- c(abc, weights_nBLP)

#    if(abs(sum(est$par[-c(1:3)]) - 1) < 0.01) break
#  }
  return(pars_nBLP)
}

get_cnBLPpars<-function(means,sds){
  weights_nBLP <- log(c(0.3, 0.3, 0.4))
  abc          <- log(c(1,1,0.1))
  pars_nBLP    <- c(abc, weights_nBLP)
  est <- optim(
    par     = pars_nBLP,
    fn      = wrapper_cnBLP,
    Y       = Y[ind_est],
    means   = means[ind_est, ],
    sds   = sds[ind_est, ],
    method  = "BFGS",
    control = list(fnscale = -1,trace = TRUE)
    
  )
  
  abc          <- exp(est$par[1:3])
  weights_nBLP <- exp(est$par[-c(1:3)])/sum(exp(est$par[-c(1:3)]))
  pars_nBLP    <- c(abc, weights_nBLP)

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
      control = list(fnscale = -1, trace = TRUE)
    )
    ab          <- exp(est$par[1:2])
    weights_CBLP <- exp(est$par[-c(1,2)])/sum(exp(est$par[-c(1,2)]))
    pars_CBLP    <- c(ab, weights_CBLP)
#    if(abs(sum(est$par[-c(1,2)]) - 1) < 0.01) break
#  }
  return(pars_CBLP)
}

get_cbcBLPpars<-function(means,sds,cab){
  weights_CBLP <- log(c(0.3, 0.3, 0.4))
  ab          <- log(c(1,1))
  pars_CBLP    <- c(ab, weights_CBLP)
  #  while(TRUE){
  
  est <- optim(
    par     = pars_CBLP,
    fn      = wrapper_cBLP,
    Y       = Y[ind_est],
    means   = means,
    sds   = sds[ind_est, ],
    cab   = cab,
    method  = "BFGS",
    control = list(fnscale = -1, trace = TRUE)
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

wrapper_SLP_nonc <- function(Y, means, sds, pars, eval = TRUE, 
                        plot = FALSE,varPIT=FALSE){
  
  con     <- pars[1]
  weights <- pars[-1]
  # normal - mean and med the same
  means_SLP <- (matrix(rep(Y, ncol(means)), ncol = ncol(means)) - means)/con
  
  ll <- log(1/con * dnorm(means_SLP, mean = 0, sd = sds) %*% weights)
  
  PIT <- pnorm(means_SLP, mean = 0, sd = sds) %*% weights
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.3), main = "SLP")
    abline( h = 1, lty = 2)
  }
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  }else{
    return(ll)
  }
}
# these two wrappers below are for plots only
wrapper_bcTLP_nonc  <- function(Y, means, sds, weights, eval = TRUE, plot = FALSE,varPIT=FALSE){
  ll <- log(dnorm(Y, mean = means, sd = sds) %*% weights)
  PIT <- pnorm(Y, mean = means, sd = sds) %*% weights
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "bcTLP")
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


wrapper_bcBLP_nonc  <- function(Y, means, sds, pars, errs, eval = TRUE, plot = FALSE,varPIT=FALSE){
  
  ab     <- pars[1:2]
  weights <- pars[-c(1,2)]
  corrected_means <- means + errs
  g      <- dnorm(Y, mean = corrected_means, sd = sds) %*% weights
  g_beta <- dbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, 
                  shape1 = ab[1], shape2 = ab[2])
  
  ll     <- log(g) + log(g_beta)
  
  PIT <- pbeta(pnorm(Y, mean = corrected_means, sd = sds) %*% weights, shape1 = ab[1], shape2 = ab[2])
  
  if(plot){
    hist(PIT, breaks = 0:10/10, col = grey(0.5), freq = FALSE, 
         ylim = c(0, 1.8), main = "bcBLP")
    abline( h = 1, lty = 2)
  } 
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}


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

wrapper_cbcBLP_nonc  <- function(Y, means, sds, pars, cab, eval = TRUE, 
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
         ylim = c(0, 1.8), main = "cbcBLP")
    abline( h = 1, lty = 2)
  } 
  if(varPIT) return(var(PIT))
  
  if(eval){
    return(sum(ll))
  } else{
    return(ll)
  }
}

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
#     control = list(fnscale = -1, trace = TRUE)
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
#     control = list(fnscale = -1, trace = TRUE)
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
#     control = list(fnscale = -1, trace = TRUE)
#   )
#   pars_indSLP <- exp(est$par)
#   return(pars_indSLP)
# }
